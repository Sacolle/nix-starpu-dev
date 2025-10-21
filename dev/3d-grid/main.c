#include <starpu.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <stdint.h>

// number of time steps that are kept in memory
#ifndef MAX_CONCURRENT_ITERATIONS
    #define MAX_CONCURRENT_ITERATIONS 100
#endif

#ifndef BASE_VOLUME_WIDTH
    #define BASE_VOLUME_WIDTH 296
#endif

#ifndef CUBE_SEGMENT_WIDTH
    #define CUBE_SEGMENT_WIDTH 30
#endif

#ifndef KERNEL_SIZE
    #define KERNEL_SIZE 4
#endif

#define VOLUME_WIDTH (BASE_VOLUME_WIDTH + KERNEL_SIZE)

#if (VOLUME_WIDTH % CUBE_SEGMENT_WIDTH) != 0
    #error A largura do volume + kernel size devem ser divisíveis pela largura do segmento.
#endif

#define WIDTH_IN_CUBES (VOLUME_WIDTH / CUBE_SEGMENT_WIDTH)

// linear idx https://stackoverflow.com/a/34363187
#define CUBE_I(x, y, z) ((x) + ((y) * CUBE_SEGMENT_WIDTH) + ((z) * CUBE_SEGMENT_WIDTH * CUBE_SEGMENT_WIDTH))
#define BLOCK_I(i, j, k) ((i) + ((j) * WIDTH_IN_CUBES) + ((k) * WIDTH_IN_CUBES * WIDTH_IN_CUBES))
#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))


#define FSTBLK(blk) ((blk) == 0)
#define LSTBLK(blk) ((blk) == (WIDTH_IN_CUBES - 1))

void print_block(double* block){
    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        printf("[\n");
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            printf("\t[");
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){
                printf("%.2f, ", block[CUBE_I(x, y, z)]);
            }
            printf("]\n");
        }
        printf("]\n");
    }
}
//initialize block, the border consists of value 0.0
void initialize_block(double* block, int i, int j, int k){
    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){

                #define START(blk, idx) (FSTBLK(blk) && (idx) < KERNEL_SIZE)
                #define END(blk, idx) (LSTBLK(blk) && (idx) >= (CUBE_SEGMENT_WIDTH - KERNEL_SIZE))
                #define INEDGE(blk_idx, cidx) (START(blk_idx, cidx) || END(blk_idx, cidx))

                //TODO: validate
                //if it starts or ends in an edge, set to 0
                if( INEDGE(k, z) || INEDGE(j, y) || INEDGE(i, x)){
                    block[CUBE_I(x, y, z)] = 0.0;
                }else{
                    block[CUBE_I(x, y, z)] = 10.0 * ((double) rand() / RAND_MAX);
                }
            }
        }
    }
}

/*
struct atomic_rc_ptr {
    //acredito que não precise de um indice interno para indexar dentro de cubes, 
    //pois isso vai ocorrer dentro do laço de execução, mas TODO: validar
    double* cubes[CUBE(WIDTH_IN_CUBES)];
    // when reaches CUBE(WIDTH_IN_CUBES)
    // then all elements in the list above where used and the list can be freed.
    // updates with the callback for the computation
    _Atomic int count; 
};

// como cada iteração tem uma forte dependência de dados na anterior
// a única forma de todos os cubos da iteração t serem usados pela iteração t + 1 é
// se todos os cubos da iteração t - 1 tiverem sido usados.
typedef struct circular_buffer_with_rc {
    struct atomic_rc_ptr  list[MAX_CONCURRENT_ITERATIONS];
    // acho q esse cara não precisa ser atômico, mas idk, n faz mal
    _Atomic int ffidx; //fist free index 
} RcCircularBuf;


RcCircularBuf rc_circular_buffer_create(){
    RcCircularBuf rcb;
    rcb.ffidx = 0;
    for(int i = 0; i < MAX_CONCURRENT_ITERATIONS; i++){
        rcb.list[i].count = 0;
        for(int j = 0; j < CUBE(WIDTH_IN_CUBES); j++){
            rcb.list[i].cubes[j] = NULL;
        }
    }
    return rcb;
}

//create buffer 
int rcb_make_buffer(RcCircularBuf* rcb, struct atomic_rc_ptr** buffer){
    // check if first free index buffer is free 
    if(rcb->list[rcb->ffidx].count != 0){
        return 1; //if not return
    }

    // if return the ptr to the buffer
    //TODO: unfinished

    rcb->ffidx = (rcb->ffidx + 1) % MAX_CONCURRENT_ITERATIONS;
}*/


typedef struct arc_data_handle {
    starpu_data_handle_t handle;
    _Atomic int uses;
} *arc_data_handle_t;

//do not call with uses < 1, please
//allocates the arc, but does not initialize the inner data handle
arc_data_handle_t create_arc(int uses){
    arc_data_handle_t ptr = (arc_data_handle_t) malloc(sizeof(struct arc_data_handle));
    if(ptr == NULL){
        return NULL;
    }
    ptr->uses = uses;
    ptr->handle = NULL;

    return ptr;
}
//if true, the underlying data can be freed safely
//fetches then subtracks, therefore if n = 1, then we n - 1, n - 1 = 0
//if true, the underlining data can be freed with no issue
bool free_arc(arc_data_handle_t ptr){
    return atomic_fetch_sub(&ptr->uses, 1) == 1;
}




struct cl_args {
    uint32_t i;
    uint32_t j;
    uint32_t k;
    uint32_t t;
};

#define MAX_NEIGHBOORS 7

//NOTE: pode reduzir o número de alocações alocando o nome e cl_args dentro de callback_args
typedef struct codelet_data {
    char name[48];
    struct cl_args args;
    arc_data_handle_t arcs[MAX_NEIGHBOORS];
} codelet_data_t;

struct codelet_data* make_codelet_data(){
    struct codelet_data* cb_args = (struct codelet_data*) malloc(sizeof(struct codelet_data));
    if(cb_args == NULL){
        return NULL;
    }
    for(int i = 0; i < MAX_NEIGHBOORS; i++){
        cb_args->arcs[i] = NULL;
    }

    return cb_args;
}

void set_cl_args(codelet_data_t* cl_data, uint32_t i, uint32_t j, uint32_t k, uint32_t t){
    cl_data->args.i = i;
    cl_data->args.j = j;
    cl_data->args.k = k;
    cl_data->args.t = t;
}

void cb_free_arcs(void* callback_arg){
    codelet_data_t* data = (codelet_data_t*) callback_arg;

    for(int i = 0; i < MAX_NEIGHBOORS; i++){
        if(data->arcs[i] != NULL){
            arc_data_handle_t arc = data->arcs[i];
            // attemps to free the data under arc, 
            // if uses == 0 after sub, it returns true and 
            //this functions needs to free the underlying data
            if(free_arc(arc)){
                starpu_data_unregister_no_coherency(arc->handle);
                free(arc);
            }
        }
    }
}

void average_filter(void *descr[], void *cl_args){

}

int main(int argc, char **argv){

    struct starpu_codelet avrg_filter_cl =
        {
            .cpu_funcs = {average_filter},
            // acredito que esse NBUFFERS vai atrapalhar na performance
            .nbuffers = STARPU_VARIABLE_NBUFFERS, 
            /*
            .modes = {
                STARPU_W, // b(x, y, z, t) resultado
                STARPU_R, // b(x, y, z, t - 1) centro
                STARPU_R, // b(x + 1, y, z, t - 1) east (leste)
                STARPU_R, // b(x - 1, y, z, t - 1) west (oeste)
                STARPU_R, // b(x, y + 1, z, t - 1) north
                STARPU_R, // b(x, y - 1, z, t - 1) south
                STARPU_R, // b(x, y, z + 1, t - 1) up
                STARPU_R, // b(x, y, z + 1, t - 1) down
            },*/
            .model = &starpu_perfmodel_nop,
            // need to set this to some value?
            // .callback_func
    };

	int ret = starpu_init(NULL);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

	assert(argc > 0);
	const int iterations = atoi(argv[1]);
    assert(iterations > 0);


    // é necessário inicializar duas matrizes iniciais, pois para computar o bloco(x,y,z,t)
    // também é necessário o bloco(x,y,z,t-1)
    // agora ao invés de um blocão, os blocos seram disjuntos
    arc_data_handle_t* prev = (arc_data_handle_t*) malloc(CUBE(WIDTH_IN_CUBES));

    // need this because the starpu automatic alocation happens on demand
    // an if can be used to free this at the end of the code or after the first iteration
    // rn it is left as a TODO: free this sometime
    double* initial_values[CUBE(WIDTH_IN_CUBES)];

    for(int k = 0; k < WIDTH_IN_CUBES; k++){
        const int neighbors_z_axis = ((k == 0) || (k == WIDTH_IN_CUBES - 1)) ? 1 : 2;

        for(int j = 0; j < WIDTH_IN_CUBES; j++){
            const int neighbors_y_axis = ((j == 0) || (j == WIDTH_IN_CUBES - 1)) ? 1 : 2;

            for(int i = 0; i < WIDTH_IN_CUBES; i++){
                const int neighbors_x_axis = ((i == 0) || (i == WIDTH_IN_CUBES - 1)) ? 1 : 2;

                const int num_of_nighbors = neighbors_x_axis + neighbors_y_axis + neighbors_z_axis;
                //cada vizinho e ele mesmo vai precisar usar
                arc_data_handle_t arc = create_arc(num_of_nighbors + 1);
                assert(arc != NULL);

                double* allocated_block = (double*) malloc(sizeof(double) * CUBE(CUBE_SEGMENT_WIDTH));
                initialize_block(allocated_block, i, j, k);

                //print_block(allocated_block);

                //let starpu allocate the data by setting home_node = -1 
                //problem: the data is only alocated when writen with STARPU_W
                //wich will only happen in the codelet, (i think)
                //so the first iterations need to allocate the memory manually
                starpu_block_data_register(&arc->handle, STARPU_MAIN_RAM, (uintptr_t) allocated_block, 
                    //stride for y      stride for z
                    CUBE_SEGMENT_WIDTH, SQUARE(CUBE_SEGMENT_WIDTH),
                    // width             height             depth
                    CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, 
                    sizeof(double)
                );
                initial_values[BLOCK_I(i, j, k)] = allocated_block;

                prev[BLOCK_I(i,j,k)] = arc; 
            }
        }
    }

    arc_data_handle_t* curr = (arc_data_handle_t*) malloc(CUBE(WIDTH_IN_CUBES));


    for(int t = 0; t < iterations; t++){
        for(int k = 0; k < WIDTH_IN_CUBES; k++){
            const int neighbors_z_axis = ((k == 0) || (k == WIDTH_IN_CUBES - 1)) ? 1 : 2;

            for(int j = 0; j < WIDTH_IN_CUBES; j++){
                const int neighbors_y_axis = ((j == 0) || (j == WIDTH_IN_CUBES - 1)) ? 1 : 2;

                for(int i = 0; i < WIDTH_IN_CUBES; i++){
                    const int neighbors_x_axis = ((i == 0) || (i == WIDTH_IN_CUBES - 1)) ? 1 : 2;

                    const int num_of_nighbors = neighbors_x_axis + neighbors_y_axis + neighbors_z_axis;
                    //cada vizinho e ele mesmo vai precisar usar
                    arc_data_handle_t arc = create_arc(num_of_nighbors + 1);
                    assert(arc != NULL);

                    //let starpu allocate the data by setting home_node = -1 
                    starpu_block_data_register(&arc->handle, -1, 0,
                        CUBE_SEGMENT_WIDTH, SQUARE(CUBE_SEGMENT_WIDTH),
                        CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, 
                        sizeof(double)
                    );

                    //add to the curr buff
                    curr[BLOCK_I(i, j, k)] = arc;

                    // 
                    struct starpu_task* task = starpu_task_create();

                    task->cl = &avrg_filter_cl;

                    struct codelet_data* cl_data = make_codelet_data();
                    assert(cl_data != NULL);

                    sprintf(cl_data->name, "[%d, %d, %d, %d]", i, j, k, t + 1);
                    task->name = cl_data->name;

                    set_cl_args(cl_data, i, j, k, t + 1);

                    task->cl_arg = &cl_data->args;
                    task->cl_arg_size = sizeof(struct cl_args);
                    // the args are alocated with the callback, so it is freed with it
                    task->cl_arg_free = 0;


                    //select the handles
                    //need to stablish an order for selection, em 3d eg dificil
                    //          ^   ^
                    //          |  /
                    //          y z
                    //          |/
                    // -- x -- >
                    // ordem ao invés de rotacional vai ser via eixo,
                    // blocos do z, depois do y, depois do x
                    // no caso primeiro vai o bloco curr, 
                    // depois o central prev e dps o resto na ordem estabelecida

                    // buffers iniciais
                    task->handles[0] = curr[BLOCK_I(i, j, k)]->handle; 
                    task->modes[0] = STARPU_W;

                    task->handles[1] = prev[BLOCK_I(i, j, k)]->handle;
                    task->modes[1] = STARPU_R;
                    //passa para o callback o arc inteiro, q a iteração anterior tem q ser liberada
                    cl_data->arcs[0] = prev[BLOCK_I(i, j, k)];

                    int nbuffers = 2;
                    
                    #define ADDBLOCK(idx) \
                        task->handles[nbuffers] = prev[idx]->handle; \
                        task->modes[nbuffers] = STARPU_R; \
                        cl_data->arcs[nbuffers - 1] = prev[idx]; \
                        nbuffers++;

                    if(!FSTBLK(k)){ // tem o bloco k - 1
                        ADDBLOCK(BLOCK_I(i, j, k - 1));
                    }
                    if(!LSTBLK(k)){ //tem o bloco k + 1
                        ADDBLOCK(BLOCK_I(i, j, k + 1));
                    }
                    if(!FSTBLK(j)){ // tem o bloco j - 1
                        ADDBLOCK(BLOCK_I(i, j - 1, k));
                    }
                    if(!LSTBLK(j)){ //tem o bloco j + 1
                        ADDBLOCK(BLOCK_I(i, j + 1, k));
                    }
                    if(!FSTBLK(i)){ // tem o bloco i - 1
                        ADDBLOCK(BLOCK_I(i - 1, j, k));
                    }
                    if(!LSTBLK(i)){ //tem o bloco i + 1
                        ADDBLOCK(BLOCK_I(i + 1, j, k));
                    }

                    assert(num_of_nighbors + 2 == nbuffers);
                    task->nbuffers = nbuffers;

                    //add the callback function
                    task->callback_func = cb_free_arcs;
                    task->callback_arg = cl_data;
                    //auto frees the args
                    task->callback_arg_free = 1;

                    ret = starpu_task_submit(task);
                    //printf("iter: %d\n, %s", ++iter,t->name);
                    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
                }
            }
        }
        prev = curr;
        // when to wait for all?
    }
    //at least after all iterations
    starpu_task_wait_for_all();

    //output the results in curr
    //desregistra os blocos e limpa as tasks
    for(int k = 0; k < WIDTH_IN_CUBES; k++){
        for(int j = 0; j < WIDTH_IN_CUBES; j++){
            for(int i = 0; i < WIDTH_IN_CUBES; i++){
                arc_data_handle_t arc = curr[BLOCK_I(i,j,k)];
                double* allocated_block = (double*) STARPU_BLOCK_GET_PTR(arc->handle);
                //print_block(allocated_block);
                
                //no need to do the ARC stuff here,
                //these blocks are not referenced by other blocks
                starpu_data_unregister(arc->handle);
                free(arc);
            }
        }
    }
	starpu_shutdown();
	return 0;
}