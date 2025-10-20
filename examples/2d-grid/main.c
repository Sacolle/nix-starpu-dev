#include <starpu.h>
#include <pthread.h>

#define IDX(i, j) ((i) + ((j)*n))
#define BLOCK(i, j) ((i) + ((j) * block_amounts_w))
#define SQ(x) ((x)*(x))

#define BELOW 0
#define RIGHT 1
#define ABOVE 2
#define LEFT 3

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

//#define DEBUG

void print_matrix(double* A, int block_amounts_w, int block_w){
    const int bs = block_w * block_w;

    #ifdef DEBUG
    for (int by = 0; by < block_amounts_w; by++){
        for (int j = 0; j < block_w; j++){
            for (int bx = 0; bx < block_amounts_w; bx++){
                for (int i = 0; i < block_w; i++){
                    const int idx = (BLOCK(bx, by) * bs) + j * block_w + i;
                    printf("%0.2f, ", A[idx]);
                }
            }
            printf("\n");
        }
    }
    printf("\n");
    #endif
}

void read_params(
    int argc, char **argv, 
    int *board_size, int *block_size, int* iterations
){
	assert(argc > 3);

	*board_size = atoi(argv[1]);
	*block_size = atoi(argv[2]);
	*iterations = atoi(argv[3]);

    assert(*iterations > 0);
}
/*
Enche a grid com os valores da borda como 0 e o centro como aleatórios
*/
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
    /*
	for (int j = 0; j < n; j++){
		for (int i = 0; i < n; i++){
            if (j < border_size || j > n - 1 - border_size || i < border_size || i > n - 1 - border_size) {
                A[IDX(i, j)] = 0.0;
                continue;
            }
			A[IDX(i, j)] = 10.0 * ((double) rand() / RAND_MAX);
		}
	}*/
}

void stencil_cpu(void *descr[], void *cl_args){
    struct starpu_variable_interface *v = (struct starpu_variable_interface *) descr[0];
    
    // um desses blocos não é passado, então segue a ordem anti-horária para o bloco omitido
	double* writeb = (double *) STARPU_VECTOR_GET_PTR(descr[0]);
	double* centerb = (double *) STARPU_VECTOR_GET_PTR(descr[1]);

    struct params *p = (struct params* ) cl_args;
    //size of the block
    const int block_amounts_w = p->block_amounts_w;
    const int n = p->block_width;
    const int block_size = SQ(n);
    const int border_size = p->border_size;

    
    const int start_x = p->block_x == 0 ? p->border_size : 0; // if this then left side
    const int end_x = n - (p->block_x == block_amounts_w - 1 ? p->border_size : 0); // if this then right side

    const int start_y = p->block_y == 0 ? p->border_size : 0; // if this then top side
    const int end_y = n - (p->block_y == block_amounts_w - 1 ? p->border_size : 0); // if this then bottom side

    // atribui o valor do bloco central ao bloco faltante
    // seguindo a conveção de primeiro para baixo, depois anti horário.
    double* adj_blocks[4];

    if(p->block_x == 0){
        // does not use left block
        if(p->block_y == 0){
            // does not use above block
            // c b3
            // b2
            adj_blocks[BELOW] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[RIGHT] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[ABOVE] = centerb;
            adj_blocks[LEFT]  = centerb;

        }else if(p->block_y == block_amounts_w - 1){
            // does not use below block
            // b3
            // c b2
            adj_blocks[BELOW] = centerb;
            adj_blocks[RIGHT] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[ABOVE] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[LEFT]  = centerb;

        }else{
            // only does not use left block
            // b4
            // c b3
            // b2
            adj_blocks[BELOW] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[RIGHT] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[ABOVE] = ( double* ) STARPU_VECTOR_GET_PTR(descr[4]);
            adj_blocks[LEFT]  = centerb;
        }
    }else if(p->block_x == block_amounts_w - 1){
        // does not use right block
        if(p->block_y == 0){
            // does not use above block
            // b3 c
            //   b2
            adj_blocks[BELOW] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[RIGHT] = centerb;
            adj_blocks[ABOVE] = centerb;
            adj_blocks[LEFT]  = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);

        }else if(p->block_y == block_amounts_w - 1){
            // does not use below block
            //   b2
            // b3 c 
            adj_blocks[BELOW] = centerb;
            adj_blocks[RIGHT] = centerb;
            adj_blocks[ABOVE] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[LEFT]  = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);

        }else{
            // only does not use right block
            //   b3
            // b4 c 
            //   b2
            adj_blocks[BELOW] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[RIGHT] = centerb;
            adj_blocks[ABOVE] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[LEFT]  = ( double* ) STARPU_VECTOR_GET_PTR(descr[4]);
        }
    }else{
        if(p->block_y == 0){
            // only does not use above block
            // b4 c b3
            //   b2
            adj_blocks[BELOW] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[RIGHT] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[ABOVE] = centerb;
            adj_blocks[LEFT]  = ( double* ) STARPU_VECTOR_GET_PTR(descr[4]);

        }else if(p->block_y == block_amounts_w - 1){
            // only does not use below block
            //   b3
            // b4 c b2
            adj_blocks[BELOW] = centerb;
            adj_blocks[RIGHT] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[ABOVE] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[LEFT]  = ( double* ) STARPU_VECTOR_GET_PTR(descr[4]);
        }else{
            //uses all blocks
            adj_blocks[BELOW] = ( double* ) STARPU_VECTOR_GET_PTR(descr[2]);
            adj_blocks[RIGHT] = ( double* ) STARPU_VECTOR_GET_PTR(descr[3]);
            adj_blocks[ABOVE] = ( double* ) STARPU_VECTOR_GET_PTR(descr[4]);
            adj_blocks[LEFT]  = ( double* ) STARPU_VECTOR_GET_PTR(descr[5]);
        }
    }
    // TODO: OS blocos estão em coordenadas locais, a conta do índice deles está errada
    for(int y = start_y; y < end_y; y++){
        for(int x = start_x; x < end_x; x++){
            
            const double below = y == n - 1 ? 
                adj_blocks[BELOW][IDX(x, 0)] : 
                centerb[IDX(x, y + 1)];

            const double above = y == 0 ? 
                adj_blocks[ABOVE][IDX(x, n - 1)] : 
                centerb[IDX(x, y - 1)];

            const double right = x == n - 1 ? 
                adj_blocks[RIGHT][IDX(0, y)] : 
                centerb[IDX(x + 1, y)];

            const double left = x == 0 ? 
                adj_blocks[LEFT][IDX(n - 1, y)] : 
                centerb[IDX(x - 1, y)];

            const int center_idx = IDX(x, y);
            writeb[center_idx] = (centerb[center_idx] + below + above + right + left) / ((double) 5.0);
        }
    }
}

struct starpu_codelet stencil6_cl =
	{
		.cpu_funcs = {stencil_cpu},
		.nbuffers = 6,
		.modes = {STARPU_W, STARPU_R, STARPU_R, STARPU_R, STARPU_R, STARPU_R},
		.model = &starpu_perfmodel_nop,
};

struct starpu_codelet stencil5_cl =
	{
		.cpu_funcs = {stencil_cpu},
		.nbuffers = 5,
		.modes = {STARPU_W, STARPU_R, STARPU_R, STARPU_R, STARPU_R},
		.model = &starpu_perfmodel_nop,
};

struct starpu_codelet stencil4_cl =
	{
		.cpu_funcs = {stencil_cpu},
		.nbuffers = 4,
		.modes = {STARPU_W, STARPU_R, STARPU_R, STARPU_R},
		.model = &starpu_perfmodel_nop,
};

void callback_func(void *callback_arg){
    pthread_t current_thread_id = pthread_self();
    printf("Callback in thread with ID %lu\n", (unsigned long)current_thread_id);
}

int main(int argc, char **argv){
	int ret = starpu_init(NULL);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

    int board_size, block_size, iterations, border_size = 1;
	read_params(argc, argv, &board_size, &block_size, &iterations);


    const int grid_size = board_size + border_size;

    assert(grid_size % block_size == 0);
    const int block_amounts_w = (grid_size / block_size);

    printf("blocks are %d\n", block_amounts_w);

	double *gridA = malloc(SQ(grid_size) * sizeof(double));
	fill(gridA, block_amounts_w, block_size, border_size);
	print_matrix(gridA, block_amounts_w, block_size);

	double *gridB = malloc(SQ(grid_size) * sizeof(double));
    memcpy(gridB, gridA, SQ(grid_size) * sizeof(double));
	print_matrix(gridB, block_amounts_w, block_size);



	starpu_data_handle_t *data_handles_A = 
        (starpu_data_handle_t*) malloc(SQ(block_amounts_w) * sizeof(*data_handles_A));

	starpu_data_handle_t *data_handles_B = 
        (starpu_data_handle_t*) malloc(SQ(block_amounts_w) * sizeof(*data_handles_B));
    
    //regitra os blocos
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
    printf("Creating the tasks...\n");
    //cria as tasks
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
                //callback
                task->starpu_task[k]->callback_func = callback_func;

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
	}

    printf("Runing the tasks over %d iterations.\n", iterations);

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

    print_matrix(gridA, block_amounts_w, block_size);
    print_matrix(gridB, block_amounts_w, block_size);

    //desregistra os blocos e limpa as tasks
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

	starpu_shutdown();
	return 0;
}
