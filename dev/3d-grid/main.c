#include <starpu.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <stdint.h>

#ifndef BASE_VOLUME_WIDTH
    #define BASE_VOLUME_WIDTH 76
#endif

#ifndef CUBE_SEGMENT_WIDTH
    #define CUBE_SEGMENT_WIDTH 8
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

struct cl_args {
    char name[48];
    uint32_t i;
    uint32_t j;
    uint32_t k;
    uint32_t t;
};

struct cl_args* make_cl_args(uint32_t i, uint32_t j, uint32_t k, uint32_t t){
    struct cl_args* cl_args = (struct cl_args*) malloc(sizeof(struct cl_args));
    if(cl_args == NULL){
        return NULL;
    }
    cl_args->i = i;
    cl_args->j = j;
    cl_args->k = k;
    cl_args->t = t;

    return cl_args;
}

void average_filter(void *descr[], void *cl_args){
    double* block_w = (double*) STARPU_BLOCK_GET_PTR( descr[0] );
    double* block_r = (double*) STARPU_BLOCK_GET_PTR( descr[1] );

    for(int i = 0; i < CUBE(CUBE_SEGMENT_WIDTH); i++){
        block_w[i] = block_r[i];
    }
}

int main(int argc, char **argv){

    struct starpu_codelet avrg_filter_cl = {
        .cpu_funcs = { average_filter },
        .nbuffers = STARPU_VARIABLE_NBUFFERS, 
        .model = &starpu_perfmodel_nop,
    };

	int ret = starpu_init(NULL);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

	assert(argc > 0);
	const int iterations = atoi(argv[1]);
    assert(iterations > 0);


    // é necessário inicializar duas matrizes iniciais, pois para computar o bloco(x,y,z,t)
    // também é necessário o bloco(x,y,z,t-1)
    // agora ao invés de um blocão, os blocos seram disjuntos
    starpu_data_handle_t* prev = (starpu_data_handle_t*) malloc(sizeof(starpu_data_handle_t) * CUBE(WIDTH_IN_CUBES));

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

                double* allocated_block = (double*) malloc(sizeof(double) * CUBE(CUBE_SEGMENT_WIDTH));
                initialize_block(allocated_block, i, j, k);

                initial_values[BLOCK_I(i, j, k)] = allocated_block;

                //print_block(allocated_block);

                //let starpu allocate the data by setting home_node = -1 
                //problem: the data is only alocated when writen with STARPU_W
                //wich will only happen in the codelet, (i think)
                //so the first iterations need to allocate the memory manually
                starpu_block_data_register(&prev[BLOCK_I(i,j,k)], STARPU_MAIN_RAM, (uintptr_t) allocated_block, 
                    //stride for y      stride for z
                    CUBE_SEGMENT_WIDTH, SQUARE(CUBE_SEGMENT_WIDTH),
                    // width             height             depth
                    CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, 
                    sizeof(double)
                );
                // assert that its properly registered
                assert(prev[BLOCK_I(i,j,k)] != NULL);
                assert(starpu_data_get_local_ptr(prev[BLOCK_I(i,j,k)]) != NULL);
            }
        }
    }

    starpu_data_handle_t* curr = (starpu_data_handle_t*) malloc(sizeof(starpu_data_handle_t) * CUBE(WIDTH_IN_CUBES));


    for(int t = 0; t < iterations; t++){
        for(int k = 0; k < WIDTH_IN_CUBES; k++){
            const int neighbors_z_axis = ((k == 0) || (k == WIDTH_IN_CUBES - 1)) ? 1 : 2;

            for(int j = 0; j < WIDTH_IN_CUBES; j++){
                const int neighbors_y_axis = ((j == 0) || (j == WIDTH_IN_CUBES - 1)) ? 1 : 2;

                for(int i = 0; i < WIDTH_IN_CUBES; i++){
                    const int neighbors_x_axis = ((i == 0) || (i == WIDTH_IN_CUBES - 1)) ? 1 : 2;

                    const int num_of_nighbors = neighbors_x_axis + neighbors_y_axis + neighbors_z_axis;
                    //cada vizinho e ele mesmo vai precisar usar

                    //add to the curr buff
                    //let starpu allocate the data by setting home_node = -1 
                    starpu_block_data_register(&curr[BLOCK_I(i, j, k)], -1, 0,
                        CUBE_SEGMENT_WIDTH, SQUARE(CUBE_SEGMENT_WIDTH),
                        CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, CUBE_SEGMENT_WIDTH, 
                        sizeof(double)
                    );

                    // 
                    struct starpu_task* task = starpu_task_create();

                    task->cl = &avrg_filter_cl;
                    
                    struct cl_args* cl_args = make_cl_args(i, j, k, t + 1);
                    assert(cl_args != NULL);

                    sprintf(cl_args->name, "[%d, %d, %d, %d]", i, j, k, t + 1);
                    task->name = cl_args->name;

                    task->cl_arg = cl_args;
                    task->cl_arg_size = sizeof(struct cl_args);
                    task->cl_arg_free = 1; // free the args after use


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
                    task->handles[0] = curr[BLOCK_I(i, j, k)]; 
                    task->modes[0] = STARPU_W;

                    task->handles[1] = prev[BLOCK_I(i, j, k)];
                    task->modes[1] = STARPU_R;

                    int nbuffers = 2;
                    
                    #define ADDBLOCK(idx) \
                        task->handles[nbuffers] = prev[idx]; \
                        task->modes[nbuffers] = STARPU_R; \
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

                    ret = starpu_task_submit(task);
                    //printf("iter: %d\n, %s", ++iter,t->name);
                    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
                }
            }
        }
        //TODO: call starpu_unregister_submit aqui para todos os valores em prev
        for(int ta = 0; ta < CUBE(WIDTH_IN_CUBES); ta++){
            starpu_data_unregister_submit(prev[ta]);
        }
        //do a swap
        starpu_data_handle_t* tmp = prev;
        prev = curr;
        curr = tmp;
        // when to wait for all?
    }
    printf("Submitted all tasks\n");
    //at least after all iterations
    starpu_task_wait_for_all();

    //TODO: explicar
    starpu_data_handle_t* result = prev;

    //output the results in curr
    //desregistra os blocos e limpa as tasks
    for(int k = 0; k < WIDTH_IN_CUBES; k++){
        for(int j = 0; j < WIDTH_IN_CUBES; j++){
            for(int i = 0; i < WIDTH_IN_CUBES; i++){
                starpu_data_handle_t handle = result[BLOCK_I(i,j,k)];

                double* result_block = initial_values[BLOCK_I(i,j,k)];
                int ret = starpu_data_peek(handle, result_block, sizeof(double) * CUBE(CUBE_SEGMENT_WIDTH));
                STARPU_CHECK_RETURN_VALUE(ret, "starpu_data_peek");

                print_block(result_block);


                
                starpu_data_unregister(handle);
            }
        }
    }
    free(prev);
    free(curr);
	starpu_shutdown();
	return 0;
}