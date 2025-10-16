#include <starpu.h>

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


struct params {
    int block_x;
    int block_y;
    int block_width;
    int block_amounts_w;
    int border_size;
};


typedef struct task_with_args {
    struct starpu_task* starpu_task[2];
    struct params params;
    char nameA[32];
    char nameB[32];
} Task;


/*
Enche a grid com os valores da borda como 0 e o centro como aleatórios
void fill(double *A, int block_amounts_w, int block_w, int border_size){
    const int bs = block_w * block_w;
    const int n = block_w;

    for (int by = 0; by < block_amounts_w; by++){
        for (int bx = 0; bx < block_amounts_w; bx++){
            for (int j = 0; j < block_w; j++){
                for (int i = 0; i < block_w; i++){
                    const int idx = (BLOCK(bx, by) * bs) + j * block_w + i;

                    // canto esquerdo
                    if(bx == 0 && i == 0){
                        A[idx] = 0.0;
                    // canto direito
                    }else if(bx == block_amounts_w - 1 && i == block_w - 1){
                        A[idx] = 0.0;
                    //topo
                    }else if(by == 0 && j == 0){
                        A[idx] = 0.0;
                    //base
                    }else if(by == block_amounts_w - 1 && j == block_w - 1){
                        A[idx] = 0.0;
                    }else{
                        A[idx] = 10.0 * ((double) rand() / RAND_MAX);
                    }

                }
            }
        }
    }
}
*/

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

    rcb->ffidx = (rcb->ffidx + 1) % MAX_CONCURRENT_ITERATIONS;
}



void stencil_cpu(void *descr[], void *cl_args){

}

int main(int argc, char **argv){

    struct starpu_codelet avrg_filter =
        {
            .cpu_funcs = {stencil_cpu},
            .nbuffers = 8, //número máximo de buffers da config default
            .modes = {
                STARPU_W, // b(x, y, z, t) resultado
                STARPU_R, // b(x, y, z, t - 1) centro
                STARPU_R, // b(x + 1, y, z, t - 1) east (leste)
                STARPU_R, // b(x - 1, y, z, t - 1) west (oeste)
                STARPU_R, // b(x, y + 1, z, t - 1) north
                STARPU_R, // b(x, y - 1, z, t - 1) south
                STARPU_R, // b(x, y, z + 1, t - 1) up
                STARPU_R, // b(x, y, z + 1, t - 1) down
            },
            .model = &starpu_perfmodel_nop,
            //add a callback?
    };

	int ret = starpu_init(NULL);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

	assert(argc > 0);
	const int iterations = atoi(argv[1]);
    assert(iterations > 0);

    //circular list to keep reference to the elements to free tem in the end
    RcCircularBuf iteration_cubes = rc_circular_buffer_create();


    /*
    // é necessário inicializar duas matrizes iniciais, pois para computar o bloco(x,y,z,t)
    // também é necessário o bloco(x,y,z,t-1)
    // agora ao invés de um blocão, os blocos seram disjuntos

	double *gridA = malloc(SQ(grid_size) * sizeof(double));
	fill(gridA, block_amounts_w, block_size, border_size);
	print_matrix(gridA, block_amounts_w, block_size);

	double *gridB = malloc(SQ(grid_size) * sizeof(double));
    memcpy(gridB, gridA, SQ(grid_size) * sizeof(double));
	print_matrix(gridB, block_amounts_w, block_size);
    */


    /*

    //data handles deverão ser lidados de forma diferente, já que agora eles são por iteração.
    
	starpu_data_handle_t *data_handles_A = 
        (starpu_data_handle_t*) malloc(SQ(block_amounts_w) * sizeof(*data_handles_A));

	starpu_data_handle_t *data_handles_B = 
        (starpu_data_handle_t*) malloc(SQ(block_amounts_w) * sizeof(*data_handles_B));
    */

    /*

    //o registro passa já para o momento de execução

    //registra os blocos
	for (int j = 0; j < block_amounts_w; j++) {
        for (int i = 0; i < block_amounts_w; i++) {
            const int block_idx = BLOCK(i, j);
            const int block_area = SQ(block_size);
            const int start_idx_of_block = block_idx * block_area;

			starpu_vector_data_register(
				&data_handles_A[block_idx],
				STARPU_MAIN_RAM,
				(uintptr_t) &(gridA[start_idx_of_block]), 
                block_area,
				sizeof(double)
			);

			starpu_vector_data_register(
				&data_handles_B[block_idx],
				STARPU_MAIN_RAM,
				(uintptr_t) &(gridB[start_idx_of_block]), 
                block_area,
				sizeof(double)
			);
		}
	}
    */
    printf("Creating the tasks...\n");
    //cria as tasks
    /*
    Task *tasks = malloc(SQ(block_amounts_w) * sizeof(Task));
    assert(tasks != NULL);

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