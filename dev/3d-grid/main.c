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


void print_block(double* block){
    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        printf("[\n");
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            printf("\t[");
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){
                printf("%.2d, ", block[CUBE_I(x, y, z)]);
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

                #define START(blk, idx) ((blk) == 0 && (idx) < KERNEL_SIZE)
                #define END(blk, idx) ((blk) == (WIDTH_IN_CUBES - 1) && (idx) >= (CUBE_SEGMENT_WIDTH - KERNEL_SIZE))
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
    arc_data_handle_t ptr = (arc_data_handle_t) malloc(sizeof(arc_data_handle_t));
    if(ptr == NULL){
        return NULL;
    }
    ptr->uses = uses;
    ptr->handle = NULL;

    return ptr;
}
//fetches then subtracks, therefore if n = 1, then we n - 1, n - 1 = 0
//if true, the underlining data can be freed with no issue
bool can_free_arc(arc_data_handle_t ptr){
    return atomic_fetch_sub(&ptr->uses, 1) == 1;
}


struct cl_args {
    uint32_t i;
    uint32_t j;
    uint32_t k;
    uint32_t t;
};

struct cl_args* make_args(uint32_t i, uint32_t j, uint32_t k, uint32_t t){
    struct cl_args* args = (struct cl_args*) malloc(sizeof(struct cl_args));
    if(args == NULL){
        return NULL;
    }
    args->i = i;
    args->j = j;
    args->k = k;
    args->t = t;

    return args;
}

//NOTE: pode reduzir o número de alocações alocando o nome e cl_args dentro de callback_args
struct callback_args {
    char* name;
    struct cl_args* args;
    arc_data_handle_t arcs[7];
};


struct callback_args* make_callback_args(char* name, struct cl_args* cl_args){
    struct callback_args* cb_args = (struct callback_args*) malloc(sizeof(struct callback_args));
    if(cb_args == NULL){
        return NULL;
    }
    cb_args->name = name;
    cb_args->args = cl_args;

    for(int i = 0; i < 7; i++){
        cb_args->arcs[i] = NULL;
    }

    return cb_args;
}


void stencil_cpu(void *descr[], void *cl_args){

}

int main(int argc, char **argv){

    struct starpu_codelet avrg_filter =
        {
            .cpu_funcs = {stencil_cpu},
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
            //add a callback?
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

                starpu_block_data_register(&arc->handle, 
                    STARPU_MAIN_RAM, -1, // -1 indica para o starpu alocar
                    //stride for y      stride for z
                    CUBE_SEGMENT_WIDTH, SQUARE(CUBE_SEGMENT_WIDTH),
                    // width             height             depth
                    CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, 
                    sizeof(double)
                );

                double* allocated_block = STARPU_BLOCK_GET_PTR(arc->handle);
                initialize_block(allocated_block, i, j, k);

                print_block(allocated_block);

                prev[BLOCK_I(i,j,k)] = arc; 
            }
        }
    }

    arc_data_handle_t* curr = (arc_data_handle_t*) malloc(CUBE(WIDTH_IN_CUBES));


    for(int t = 1; t < iterations; t++){
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

                    starpu_block_data_register(&arc->handle, 
                        STARPU_MAIN_RAM, -1,
                        CUBE_SEGMENT_WIDTH, SQUARE(CUBE_SEGMENT_WIDTH),
                        CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, 
                        sizeof(double)
                    );

                    // 
                    struct starpu_task* task = starpu_task_create();

                    task->cl = &stencil_cpu;

                    char* name = (char*) malloc(64);
                    assert(name != NULL);
                    sprintf(name, "[%d, %d, %d, %d]", i, j, k, t);

                    task->name = name;

                    struct cl_args* args = make_args(i,j,k,t);
                    assert(args != NULL);

                    task->cl_arg = args;
                    task->cl_arg_size = sizeof(struct cl_args);

                    struct callback_args* cb_agrs = make_callback_args(name, args);
                    assert(cb_agrs != NULL);

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
                    task->handles[0] = curr[BLOCK_I(i, j, k)]->handle; //TODO: allocate the curr
                    task->modes[0] = STARPU_W;

                    task->handles[1] = prev[BLOCK_I(i, j, k)]->handle;
                    task->modes[1] = STARPU_R;
                    //passa para o callback o arc inteiro, q a iteração anterior tem q ser liberada
                    cb_agrs->arcs[0] = prev[BLOCK_I(i, j, k)];

                    int nbuffers = 2;
                    
                    if(k != 0){
                        // tem o bloco k - 1
                        const int idx = BLOCK_I(i, j, k - 1);
                        task->handles[nbuffers] = prev[idx]->handle;
                        task->modes[nbuffers] = STARPU_R;
                        cb_agrs->arcs[nbuffers - 1] = prev[idx];

                        nbuffers++;
                    }

                    if(k != WIDTH_IN_CUBES - 1){
                        //tem o bloco k + 1
                        const int idx = BLOCK_I(i, j, k + 1);
                        task->handles[nbuffers] = prev[idx]->handle;
                        task->modes[nbuffers] = STARPU_R;
                        cb_agrs->arcs[nbuffers - 1] = prev[idx];

                        nbuffers++;
                    }

                    if(j != 0){
                        // tem o bloco j - 1
                        const int idx = BLOCK_I(i, j - 1, k);
                        task->handles[nbuffers] = prev[idx]->handle;
                        task->modes[nbuffers] = STARPU_R;
                        cb_agrs->arcs[nbuffers - 1] = prev[idx];

                        nbuffers++;
                    }

                    if(j != WIDTH_IN_CUBES - 1){
                        //tem o bloco j + 1
                        const int idx = BLOCK_I(i, j + 1, k);
                        task->handles[nbuffers] = prev[idx]->handle;
                        task->modes[nbuffers] = STARPU_R;
                        cb_agrs->arcs[nbuffers - 1] = prev[idx];

                        nbuffers++;
                    }

                    if(i != 0){
                        // tem o bloco j - 1
                        const int idx = BLOCK_I(i - 1, j, k);
                        task->handles[nbuffers] = prev[idx]->handle;
                        task->modes[nbuffers] = STARPU_R;
                        cb_agrs->arcs[nbuffers - 1] = prev[idx];

                        nbuffers++;
                    }

                    if(i != WIDTH_IN_CUBES - 1){
                        //tem o bloco j + 1
                        const int idx = BLOCK_I(i + 1, j, k);
                        task->handles[nbuffers] = prev[idx]->handle;
                        task->modes[nbuffers] = STARPU_R;
                        cb_agrs->arcs[nbuffers - 1] = prev[idx];

                        nbuffers++;
                    }

                    assert(num_of_nighbors + 2 == nbuffers);
                    task->nbuffers = nbuffers;

                    //TODO: submit
                }
            }
        }
        prev = curr;
    }
    /*
	for (int j = 0; j < block_amounts_w; j++) {
        for (int i = 0; i < block_amounts_w; i++) {
            Task* task = &tasks[BLOCK(i,j)];

            task->starpu_task[0] = starpu_task_create();
            task->starpu_task[1] = starpu_task_create();

            sprintf(task->nameA, "b(%d, %d)A", i, j);
            task->starpu_task[0]->name = task->nameA;

            sprintf(task->nameB, "b(%d, %d)B", i, j);
            task->starpu_task[1]->name = task->nameB;

            task->params.block_width = block_size;
            task->params.block_amounts_w = block_amounts_w;
            task->params.border_size = border_size;
            task->params.block_x = i;
            task->params.block_y = j;

            // Ordem anti-horária
            //    b4
            // b5 b1 b3  ->  b0
            //    b2
            for(int k = 0; k < 2; k++){
                task->starpu_task[k]->cl_arg = &task->params;
                task->starpu_task[k]->cl_arg_size = sizeof(struct params);

                // assing the write and the read buffers
                //write from one data_handle to the other
                if(k == 0){
                    // write to B while reading from A
                    task->starpu_task[0]->handles[0] = data_handles_B[BLOCK(i, j)];
                    task->starpu_task[0]->modes[0] = STARPU_W;

                    task->starpu_task[0]->handles[1] = data_handles_A[BLOCK(i, j)];
                    task->starpu_task[0]->modes[1] = STARPU_R;
                }else{
                    // write to A while reading from B
                    task->starpu_task[1]->handles[0] = data_handles_A[BLOCK(i, j)];
                    task->starpu_task[1]->modes[0] = STARPU_W;

                    task->starpu_task[1]->handles[1] = data_handles_B[BLOCK(i, j)];
                    task->starpu_task[1]->modes[1] = STARPU_R;
                }

                int nbuffers = 2;

                starpu_data_handle_t *data_handle = k == 0 ? data_handles_A : data_handles_B;
                //below
                if(j != block_amounts_w - 1){
                    task->starpu_task[k]->handles[nbuffers] = data_handle[BLOCK(i, j + 1)];
                    task->starpu_task[k]->modes[nbuffers] = STARPU_R;
                    nbuffers++;
                }
                //right
                if(i != block_amounts_w - 1){
                    task->starpu_task[k]->handles[nbuffers] = data_handle[BLOCK(i + 1, j)];
                    task->starpu_task[k]->modes[nbuffers] = STARPU_R;
                    nbuffers++;
                }
                //above
                if(j != 0){
                    task->starpu_task[k]->handles[nbuffers] = data_handle[BLOCK(i, j - 1)];
                    task->starpu_task[k]->modes[nbuffers] = STARPU_R;
                    nbuffers++;
                }
                //left
                if(i != 0){
                    task->starpu_task[k]->handles[nbuffers] = data_handle[BLOCK(i - 1, j)];
                    task->starpu_task[k]->modes[nbuffers] = STARPU_R;
                    nbuffers++;
                } 
                task->starpu_task[k]->nbuffers = nbuffers;

                switch (nbuffers) {
                case 6:
                    task->starpu_task[k]->cl = &stencil6_cl; 
                    break;
                case 5:
                    task->starpu_task[k]->cl = &stencil5_cl; 
                    break;
                case 4:
                    task->starpu_task[k]->cl = &stencil4_cl; 
                    break;
                
                default:
                    printf("ERRO: numero de buffers errado");
                    exit(1);
                    break;
                }
            }
		}
	}*/

    printf("Runing the tasks over %d iterations.\n", iterations);
    /*
    int iter = 0;
    while(iterations--){
        for (int j = 0; j < block_amounts_w; j++) {
            for (int i = 0; i < block_amounts_w; i++) {
                const int task_idx = BLOCK(i, j);
                // TODO: checar se na duplicação ele não libera a task automaticamente no fim
                // aka isso fica com um puta vazamento de memória.
                struct starpu_task* t;
                if(iterations % 2 == 0){
                    t = starpu_task_dup(tasks[task_idx].starpu_task[0]);
                    //printf("Running task read A write B\n");
                }else{
                    t = starpu_task_dup(tasks[task_idx].starpu_task[1]);
                    //printf("Running task read B write A\n");
                }
                ret = starpu_task_submit(t);
                //printf("iter: %d\n, %s", ++iter,t->name);
                STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
            }
        }
        starpu_task_wait_for_all();
    }
    */

    //desregistra os blocos e limpa as tasks
    /*
	for (int j = 0; j < block_amounts_w; j++) {
        for (int i = 0; i < block_amounts_w; i++) {
            const int block = BLOCK(i, j);

            starpu_task_clean((tasks[block]).starpu_task[0]);
            free((tasks[block]).starpu_task[0]);

            starpu_task_clean((tasks[block]).starpu_task[1]);
            free((tasks[block]).starpu_task[1]);

			starpu_data_unregister_submit(data_handles_A[block]);
			starpu_data_unregister_submit(data_handles_B[block]);
		}
	}
    free(tasks);
    */
	starpu_shutdown();
	return 0;
}