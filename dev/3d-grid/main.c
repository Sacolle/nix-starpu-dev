#include <starpu.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <stdint.h>

const int BASE_VOLUME_WIDTH = 5000;
const int CUBE_SEGMENT_WIDTH = 139;
const int KERNEL_SIZE = 4;
const int VOLUME_WIDTH = (BASE_VOLUME_WIDTH + KERNEL_SIZE);
const int WIDTH_IN_CUBES = (VOLUME_WIDTH / CUBE_SEGMENT_WIDTH);

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

#define START(blk, idx) (FSTBLK(blk) && (idx) < KERNEL_SIZE)
#define END(blk, idx) (LSTBLK(blk) && (idx) >= (CUBE_SEGMENT_WIDTH - KERNEL_SIZE))
#define INEDGE(blk_idx, cidx) (START(blk_idx, cidx) || END(blk_idx, cidx))
//initialize block, the border consists of value 0.0
void initialize_block(double* block, int i, int j, int k){
    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){
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

void clear_block(double* block){
    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){
                    block[CUBE_I(x, y, z)] = 0.0;
            }
        }
    }
}

void clear_pointers(starpu_data_handle_t* list){
    for(int i = 0; i < CUBE(WIDTH_IN_CUBES); i++){
        list[i] = NULL;
    }
}

#define EPSILON 0.0001

void assert_block_clear_edge(double* block, int i, int j, int k){
    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){
                const int idx = CUBE_I(x, y, z);
                if( INEDGE(k, z) || INEDGE(j, y) || INEDGE(i, x)){
                    assert(block[idx] < EPSILON);
                    assert(block[idx] > -EPSILON);
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

    struct cl_args* args = (struct cl_args*) cl_args;
    const uint32_t i = args->i;
    const uint32_t j = args->j;
    const uint32_t k = args->k;

    double *const block_w = (double*) STARPU_BLOCK_GET_PTR( descr[0] );
    double *const block_r = (double*) STARPU_BLOCK_GET_PTR( descr[1] );

    double *const z_minus = (double*) STARPU_BLOCK_GET_PTR( descr[2] );
    double *const z_plus = (double*) STARPU_BLOCK_GET_PTR( descr[3] );

    double *const y_minus = (double*) STARPU_BLOCK_GET_PTR( descr[4] );
    double *const y_plus = (double*) STARPU_BLOCK_GET_PTR( descr[5] );

    double *const x_minus = (double*) STARPU_BLOCK_GET_PTR( descr[6] );
    double *const x_plus = (double*) STARPU_BLOCK_GET_PTR( descr[7] );

    for(int z = 0; z < CUBE_SEGMENT_WIDTH; z++){
        for(int y = 0; y < CUBE_SEGMENT_WIDTH; y++){
            for(int x = 0; x < CUBE_SEGMENT_WIDTH; x++){
                // need to set the border to 0, because the cube might contain tash data
                if( INEDGE(k, z) || INEDGE(j, y) || INEDGE(i, x)){
                    block_w[CUBE_I(x,y,z)] = 0.0;
                }else {
                    double res = 0.0;
                    for(int k = -KERNEL_SIZE; k <= KERNEL_SIZE; k++){ 
                        const int cell_idx = z + k; 
                        
                        if(cell_idx < 0){ 
                            res += z_minus[CUBE_I(x, y, CUBE_SEGMENT_WIDTH + cell_idx)]; 
                        }else if(cell_idx >= CUBE_SEGMENT_WIDTH){ 
                            res += z_plus[CUBE_I(x, y, cell_idx - CUBE_SEGMENT_WIDTH)]; 
                        }else{ 
                            res += block_r[CUBE_I(x, y, cell_idx)]; 
                        } 
                    }

                    for(int k = -KERNEL_SIZE; k <= KERNEL_SIZE; k++){ 
                        const int cell_idx = y + k; 
                        if(cell_idx < 0){ 
                            res += y_minus[CUBE_I(x, CUBE_SEGMENT_WIDTH + cell_idx, z)]; 
                        }else if(cell_idx >= CUBE_SEGMENT_WIDTH){ 
                            res += y_plus[CUBE_I(x, cell_idx - CUBE_SEGMENT_WIDTH, z)]; 
                        }else{ 
                            res += block_r[CUBE_I(x, cell_idx, z)]; 
                        } 
                    }

                    for(int k = -KERNEL_SIZE; k <= KERNEL_SIZE; k++){ 
                        const int cell_idx = x + k; 
                        if(cell_idx < 0){ 
                            res += x_minus[CUBE_I(CUBE_SEGMENT_WIDTH + cell_idx, y, z)]; 
                        }else if(cell_idx >= CUBE_SEGMENT_WIDTH){ 
                            res += x_plus[CUBE_I(cell_idx - CUBE_SEGMENT_WIDTH, y, z)]; 
                        }else{ 
                            res += block_r[CUBE_I(cell_idx, y, z)]; 
                        } 
                    }
                        
                    block_w[CUBE_I(x, y, z)] = (res - (2 * block_r[CUBE_I(x, y, z)])) / 25.0;
                }
            }
        }
    }
    assert_block_clear_edge(block_w, i, j, k);

}

int main(int argc, char **argv){
  if (VOLUME_WIDTH % CUBE_SEGMENT_WIDTH != 0) {
    fprintf(stderr, "A largura do volume + kernel size devem ser divisíveis pela largura do segmento.\n");
    return 1;
  }


  
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
    clear_pointers(prev);

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

                //print_block(allocated_block);

                initial_values[BLOCK_I(i, j, k)] = allocated_block;

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
    clear_pointers(curr);

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
                    sprintf(cl_args->name, "avg", i, j, k, t + 1);
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
                    // blocos do z (-1, +1), depois do y, depois do x
                    // no caso primeiro vai o bloco curr, 
                    // depois o central prev e dps o resto na ordem estabelecida

                    // buffers iniciais
                    task->handles[0] = curr[BLOCK_I(i, j, k)]; 
                    task->modes[0] = STARPU_W;

                    task->handles[1] = prev[BLOCK_I(i, j, k)];
                    task->modes[1] = STARPU_R;

                    int nbuffers = 2;

                    // else add at that position the block itself
                    #define ADDBLOCKIF(cond, idx) \
                        task->handles[nbuffers] = (cond) ? prev[idx] : prev[BLOCK_I(i,j,k)]; \
                        task->modes[nbuffers] = STARPU_R; \
                        nbuffers++; 
                    
                    ADDBLOCKIF(!FSTBLK(k), BLOCK_I(i, j, k - 1));
                    ADDBLOCKIF(!LSTBLK(k), BLOCK_I(i, j, k + 1));
                    ADDBLOCKIF(!FSTBLK(j), BLOCK_I(i, j - 1, k));
                    ADDBLOCKIF(!LSTBLK(j), BLOCK_I(i, j + 1, k));
                    ADDBLOCKIF(!FSTBLK(i), BLOCK_I(i - 1, j, k));
                    ADDBLOCKIF(!LSTBLK(i), BLOCK_I(i + 1, j, k));

                    //assert(num_of_nighbors + 2 == nbuffers);
                    assert(nbuffers == 8);
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
        //write on a clear curr
        clear_pointers(curr);
        // when to wait for all?
    }
    printf("Submitted all tasks\n");
    //at least after all iterations
    starpu_task_wait_for_all();

    //TODO: explicar
    starpu_data_handle_t* result = prev;

    //output the results in curr
    //desregistra os blocos e limpa as tasks
    double* result_block = (double*) malloc(sizeof(double) * CUBE(CUBE_SEGMENT_WIDTH));
    clear_block(result_block);

    for(int k = 0; k < WIDTH_IN_CUBES; k++){
        for(int j = 0; j < WIDTH_IN_CUBES; j++){
            for(int i = 0; i < WIDTH_IN_CUBES; i++){
                starpu_data_handle_t handle = result[BLOCK_I(i,j,k)];

                // first need to acquire the data
                starpu_data_acquire(handle, STARPU_R);
                //double* result_block = initial_values[BLOCK_I(i,j,k)];
                starpu_ssize_t og_size = sizeof(double) * CUBE(CUBE_SEGMENT_WIDTH);
                starpu_ssize_t size = og_size;
                int ret = starpu_data_pack(handle, (void**)&result_block, &size);
                assert(og_size == size);
                STARPU_CHECK_RETURN_VALUE(ret, "starpu_data_peek");

                print_block(result_block);

               //assert_block_clear_edge(result_block, i, j, k);

                clear_block(result_block);
                //then release
                starpu_data_release(handle);
                starpu_data_unregister(handle);
            }
        }
    }
    free(prev);
    free(curr);
	starpu_shutdown();
	return 0;
}