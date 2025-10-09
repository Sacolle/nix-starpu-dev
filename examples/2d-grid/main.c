#include <starpu.h>
#define IDX(i, j) ((i) + ((j)*n))
#define BLOCK(i, j) ((i) + ((j) * block_amounts_w))
#define SQ(x) ((x)*(x))

struct params {
    int block_x;
    int block_y;
    int block_width;
    int block_amounts_w;
    int border_size;
};

typedef struct task_with_args {
    struct starpu_task* starpu_task;
    struct params params;
    char name[32];
} Task;

void print_matrix(double* A, int n){
	for (int j = 0; j < n; j++){
		for (int i = 0; i < n; i++){
			printf("%0.1f, ", A[IDX(i, j)]);
		}
		printf("\n");
	}
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
void fill(double *A, int n, int border_size){
	for (int j = 0; j < n; j++){
		for (int i = 0; i < n; i++){
            if (j < border_size || j > n - 1 - border_size || i < border_size || i > n - 1 - border_size) {
                A[IDX(i, j)] = 0.0;
                continue;
            }
			A[IDX(i, j)] = 10.0 * ((double) rand() / RAND_MAX);
		}
	}
}

void stencil5_cpu(void *descr[], void *cl_args){
    struct params *p = (struct params* ) cl_args;

    /*                      
	double* b0 = (double *) STARPU_VARIABLE_GET_PTR(descr[0]);
	double* b1 = (double *) STARPU_VARIABLE_GET_PTR(descr[1]);
	double* b2 = (double *) STARPU_VARIABLE_GET_PTR(descr[2]);
	double* b3 = (double *) STARPU_VARIABLE_GET_PTR(descr[3]);
	double* b4 = (double *) STARPU_VARIABLE_GET_PTR(descr[4]);

    //size of the block
    const int block_amounts_w = p->block_amounts_w;
    const int n = p->block_width;

    //const int block_i = p->block_x + p->block_y * p->block_width;
    //const int block_start = block_size * block_i;
    
    const int start_x = p->block_x == 0 ? p->border_size : 0;
    const int end_x = n - (p->block_x == block_amounts_w - 1 ? p->border_size : 0);

    const int start_y = p->block_y == 0 ? p->border_size : 0;
    const int end_y = n - (p->block_y == block_amounts_w - 1 ? p->border_size : 0);


    for(int y = start_y; y < end_y; y++){
        for(int x = start_x; x < end_x; x++){
            
            const double below = y == n - 1 ? 
                b1[BLOCK(p->block_x, p->block_y + 1) * SQ(n) + IDX(x, 0)] : 
                b0[BLOCK(p->block_x, p->block_y) * SQ(n) + IDX(x, y + 1)];
            // n - 1 é a última linha do bloco, que é o overflow 
            // quando se está voltando uma linha
            // a apartir da primeira linha do bloco abaixo
            const double above = y == 0 ? 
                b3[BLOCK(p->block_x, p->block_y - 1) * SQ(n) + IDX(x, n - 1)] : 
                b0[BLOCK(p->block_x, p->block_y) * SQ(n) + IDX(x, y - 1)];

            const double right = x == n - 1 ? 
                b2[BLOCK(p->block_x + 1, p->block_y) * SQ(n) + IDX(0, y)] : 
                b0[BLOCK(p->block_x, p->block_y) * SQ(n) + IDX(x + 1, y)];

            const double left = x == 0 ? 
                b4[BLOCK(p->block_x - 1, p->block_y) * SQ(n) + IDX(n - 1, y)] : 
                b0[BLOCK(p->block_x, p->block_y) * SQ(n) + IDX(x - 1, y)];

            const int center_idx = BLOCK(p->block_x, p->block_y) * SQ(n) + IDX(x, y);
            b0[center_idx] = (b0[center_idx] + below + above + right + left) / 5.0;
        }
    }
    */
}


struct starpu_codelet stencil5_cl =
	{
		.cpu_funcs = {stencil5_cpu},
		.nbuffers = 5,
		.modes = {STARPU_RW, STARPU_R, STARPU_R, STARPU_R, STARPU_R},
		.model = &starpu_perfmodel_nop,
};

struct starpu_codelet stencil4_cl =
	{
		.cpu_funcs = {stencil5_cpu},
		.nbuffers = 4,
		.modes = {STARPU_RW, STARPU_R, STARPU_R, STARPU_R},
		.model = &starpu_perfmodel_nop,
};

struct starpu_codelet stencil3_cl =
	{
		.cpu_funcs = {stencil5_cpu},
		.nbuffers = 3,
		.modes = {STARPU_RW, STARPU_R, STARPU_R},
		.model = &starpu_perfmodel_nop,
};

int main(int argc, char **argv){
	int ret = starpu_init(NULL);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

    int board_size, block_size, iterations, border_size = 1;
	read_params(argc, argv, &board_size, &block_size, &iterations);

    const int grid_size = board_size + border_size;
	double *grid = malloc(SQ(grid_size) * sizeof(double));
	fill(grid, grid_size, border_size);

	print_matrix(grid, grid_size);
    printf("\n");

    assert(grid_size % block_size == 0);
    const int block_amounts_w = (grid_size / block_size);

    printf("blocks are %d\n", block_amounts_w);

	starpu_data_handle_t *data_handles = 
        (starpu_data_handle_t*) malloc(SQ(block_amounts_w) * sizeof(*data_handles));
    
    //regitra os blocos
	for (int j = 0; j < block_amounts_w; j++) {
        for (int i = 0; i < block_amounts_w; i++) {
            const int block_idx = BLOCK(i, j);
            const int block_area = SQ(block_size);
            const int start_idx_of_block = block_idx * block_area;

			starpu_vector_data_register(
				&data_handles[block_idx],
				STARPU_MAIN_RAM,
				(uintptr_t) &(grid[start_idx_of_block]), 
                block_area,
				sizeof(double)
			);
		}
	}
    //cria as tasks
    Task *tasks = malloc(SQ(block_amounts_w) * sizeof(Task));
    assert(tasks != NULL);

	for (int j = 0; j < block_amounts_w; j++) {
        for (int i = 0; i < block_amounts_w; i++) {
            Task* task = &tasks[BLOCK(i,j)];
            task->starpu_task = malloc(sizeof(struct starpu_task));
            assert(task->starpu_task != NULL);
            starpu_task_init(task->starpu_task);

            sprintf(task->name, "bloco(%d, %d)", i, j);
            task->starpu_task->name = task->name;

            task->params.block_width = block_size;
            task->params.block_amounts_w = block_amounts_w;
            task->params.border_size = border_size;
            task->params.block_x = i;
            task->params.block_y = j;

            task->starpu_task->cl_arg = &task->params;
            task->starpu_task->cl_arg_size = sizeof(struct params);



            // Ordem anti-horária
            //    b3
            // b4 b0 b2
            //    b1
            //b0
            task->starpu_task->handles[0] = data_handles[BLOCK(i, j)];
            task->starpu_task->modes[0] = STARPU_RW;

            int nbuffers = 1;
            //b1
            if(j != block_amounts_w - 1){
                task->starpu_task->handles[nbuffers] = data_handles[BLOCK(i, j + 1)];
                task->starpu_task->modes[nbuffers] = STARPU_R;
                nbuffers++;
            }
            //b2
            if(i != block_amounts_w - 1){
                task->starpu_task->handles[nbuffers] = data_handles[BLOCK(i + 1, j)];
                task->starpu_task->modes[nbuffers] = STARPU_R;
                nbuffers++;
            }
            //b3
            if(j != 0){
                task->starpu_task->handles[nbuffers] = data_handles[BLOCK(i, j - 1)];
                task->starpu_task->modes[nbuffers] = STARPU_R;
                nbuffers++;
            }
            //b4
            if(i != 0){
                task->starpu_task->handles[nbuffers] = data_handles[BLOCK(i - 1, j)];
                task->starpu_task->modes[nbuffers] = STARPU_R;
                nbuffers++;
            } 
            task->starpu_task->nbuffers = nbuffers;

            switch (nbuffers) {
            case 5:
                task->starpu_task->cl = &stencil5_cl; 
                break;
            case 4:
                task->starpu_task->cl = &stencil4_cl; 
                break;
            case 3:
                task->starpu_task->cl = &stencil3_cl; 
                break;
            
            default:
                printf("ERRO: numero de buffers errado");
                exit(1);
                break;
            }
		}
	}

    int iter = 0;
    while(iterations--){
        for (int j = 0; j < block_amounts_w; j++) {
            for (int i = 0; i < block_amounts_w; i++) {
                const int task_idx = BLOCK(i, j);
                ret = starpu_task_submit(tasks[task_idx].starpu_task);
                printf("iter: %d\n", ++iter);
                STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
            }
        }
        starpu_task_wait_for_all();
    }

    //desregistra os blocos e limpa as tasks
	for (int j = 0; j < block_amounts_w; j++) {
        for (int i = 0; i < block_amounts_w; i++) {
            const int block = BLOCK(i, j);

            starpu_task_clean((tasks[block]).starpu_task);
            free((tasks[block]).starpu_task);

			starpu_data_unregister(data_handles[block]);
		}
	}
    free(tasks);
	print_matrix(grid, grid_size);

	starpu_shutdown();
	return 0;
}