#include <mpi.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <cstdlib> // for rand() and srand()
#include <ctime> // for time()
#define NUM_CITIES 194
#define POPULATION_SIZE 200
#define MUTATION_RATIO 40
#define NUM_GENERATION 5000
#define MIGRATION_PERIOD 100

class City{
public:
    int x;
    int y;
    void print() {
        std::cout <<"(" << x << "," << y << ")";
    }
};

class Country{
public:
    int num_of_cities;
    City *city_p;
};

class Path{
public:
    float fitness;
    float distance;
    int citylist[NUM_CITIES];
    int i;
    float probability;
    void print() {
        for (i = 0; i < NUM_CITIES; i++) {
        std::cout << citylist[i]<<"-";
        }
    }
};

class Population{
public:
    Path *pa;
    int no_of_path;
    float total_fit;
};

float create_random(int seed){
    srand(seed*time(0)); 
    float random_number = rand() / (float)RAND_MAX;
    return random_number;
}

int load_city_text(char **file){

    FILE * text_file;char *data;long file_size;
        
    text_file = fopen("qa194.tsp.txt", "rb");
    // get the file size
    fseek (text_file , 0 , SEEK_END);
    file_size = ftell (text_file);
    rewind (text_file);
    // malloc help to get memory to load the file:
    data = (char*) malloc (sizeof(char)*file_size);    
    // read the file to data:
    fread (data,1,file_size,text_file);
    *file = data;
    fclose (text_file);   
}
int feed_data(Country *coun){

    char *t_file = NULL;char *word = NULL;char *position = NULL;
    int i, j;

    load_city_text(&t_file);
    position = strstr(t_file, "NODE_COORD_SECTION");
    position += 19;
    coun->num_of_cities = NUM_CITIES;
    coun->city_p= (City*)malloc(sizeof(City)*coun->num_of_cities);
    word = position;
    i = 0;
    while(word != NULL && i < NUM_CITIES){
        // don't store the city name, as the city name at position 0
        if(i==0)
            word = strtok(position, " ");
        else
            word = strtok(NULL, " ");

        word = strtok(NULL, " ");
        coun->city_p[i].x = atoi(word);
        word = strtok(NULL, " \n");
        coun->city_p[i].y = atoi(word);
        if(*word == '\n'){
            word++;
        }
        i++;
    }
    free(t_file);
}

float cal_Fit(int *city_l, Country *coun){
    // Input: the first city of a path and country
    // Output: the fitness of this path
    float distance = 0.0;
    float fitness;
    int i, City1, City2;
    int xd, yd;
    for (i=0; i< NUM_CITIES -1; i++){
        City1 = city_l[i];
        City2 = city_l[i+1];        
        distance = distance + (sqrt(pow((coun->city_p[City1].x - coun->city_p[City2].x),2) +pow((coun->city_p[City1].y - coun->city_p[City2].y),2)));       
    }
    // Need to add the distance from the "end city" to the " start city" to make a circle
    City1 = city_l[0];
    City2 = city_l[NUM_CITIES -1]; 
    distance = distance + (sqrt(pow((coun->city_p[City1].x - coun->city_p[City2].x),2) +pow((coun->city_p[City1].y - coun->city_p[City2].y),2)));
    return distance;
}

int create_Path(Path *path_p, int seed){
    // First, create a temporary list of cities
    int j;
    for(j=0; j<NUM_CITIES; j++){
        path_p->citylist[j] = j; 
    }
    // Then randomly select a position, then swap it with the current position.
    int i;    
    for (i = NUM_CITIES-1; i > 0; i--) {
        int j = (unsigned int) ((create_random(seed)+drand48())*(i+1)/2);
        int t = path_p->citylist[j];
        path_p->citylist[j] = path_p->citylist[i];
        path_p->citylist[i] = t;
    }
}

int create_Pop(Population *popu, int num, int seed){
    popu->no_of_path = num;
    popu->pa = (Path*)malloc(num * sizeof(Path));
    int i;
    for(i=0; i<num; i++){
        create_Path(&popu->pa[i],seed);
        }
}

int cal_fit_Pop(Population *popul, Country *count){
    int i;
    float total_Fit = 0.0;
    for(i=0; i<popul->no_of_path; i++){
        popul->pa[i].distance = cal_Fit(&popul->pa[i].citylist[0], count);
        // The return type is int so multiple by a big number help when the distance differ by a small fraction after dot
        popul->pa[i].fitness = 1000000/popul->pa[i].distance;
        total_Fit = total_Fit + popul->pa[i].fitness;
    }
    popul->total_fit = total_Fit;
}

int cmpfunc (const void *a, const void *b){
  Path * P1 = (Path *)a;
  Path * P2 = (Path *)b;  
  return (P1->distance - P2->distance);
}

float sort_fitness(Population *popu){
    qsort(popu->pa, popu->no_of_path, sizeof(Path), cmpfunc);
}

void update_probability(Population *popu){
    int i;
    for (i=0; i < popu->no_of_path; i++){
        popu->pa[i].probability = popu->pa[i].fitness/popu->total_fit;
    }
}

void print_pop(Population *popu){
    printf("The current best distance is:");
    //popu->pa[0].print();
    printf("\n");
    printf("%f \n",popu->pa[0].distance);
    printf("-------------------- \n");
}

float Roulette_Wheel_Selection(Population *popu){
    float sum;
    sum = 0;
    float r = drand48();
    int i;
    for (i=0; i < popu->no_of_path; i++)
        {
            sum = sum + popu->pa[i].probability;
            if (r < sum) {
                return i;
            }
        }
}

float inherite_child(int OX,int *child,Path *parent1,Path *parent2){
    int i;
    for (i = 0; i < OX; i++ ){
        child[i] = parent1->citylist[i];
    }
    //The remain city in the son need to fill up after copy some parts from father
    int remain = NUM_CITIES - OX; 
    int count = 0; 
    for (int i=0; i<NUM_CITIES; i++) {
        bool exist = false;
        for (int j=0; j<OX; j++) {
            // If the city is in the child, exit this loop
            if (child[j] == parent2->citylist[i]){
                exist = true;
                break;
            }
        }            
        // If the city is not exist in the child, then copy from mother to child
        if (!exist){         
            child[OX+count] = parent2->citylist[i];
            count = count + 1;
        }
            
        // Stop when fill up all the cities
        if (count == remain) break;
    }
}

// Doing one point cross over. Randomly select a cross over point
// Then copy all from 0 to that point of one parent, and the remaining from another parent
// with condition that the remaining is not in the child. 
float order_crossover(Population *popu){
    int mother;
    int father;
    mother = Roulette_Wheel_Selection(popu);
    father = Roulette_Wheel_Selection(popu);
    while (father == mother){
        father = Roulette_Wheel_Selection(popu);
    }
    int OX ;
    int temp;
    int i;
    // Generate random number OX in range 1 to NUM_CITIES - 1, run the random several times.
    for (i=0; i < 3; i++){
        OX = (NUM_CITIES -1)*(double)rand()/RAND_MAX + 1;
    }
    
    int son[NUM_CITIES] = {0};
    int daughter[NUM_CITIES] = {0};
    inherite_child(OX,son,&popu->pa[father],&popu->pa[mother]);
    inherite_child(OX,daughter,&popu->pa[mother],&popu->pa[father]);
  
    for (i=0; i<NUM_CITIES; i++){
            popu->pa[popu->no_of_path-1].citylist[i] = son[i];
            popu->pa[popu->no_of_path-2].citylist[i] = daughter[i];
        }
}

float mutation(Population *popu){
    // rand()%100 create a random number between 0 to 99
    if(rand()%100 < MUTATION_RATIO){

        int r;
        // r is a random number between 2 and population size - 1, because we don't want to mutate the 2 best paths
        r = 2+rand() % (popu->no_of_path-2) ;
        int i; 
        for (i = NUM_CITIES-1; i > 0; i--) {
            int j = (unsigned int) (drand48()*(i+1));
            int t = popu->pa[r].citylist[j];
            popu->pa[r].citylist[j] = popu->pa[r].citylist[i];
            popu->pa[r].citylist[i] = t;
        }
    }
}
int main(int argc, char **argv)
{
    
    int rank; int size; 
    int i; int j; float dist ; float fit;
    srand (time(NULL));
    double start_time, end_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);  
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    Country CT;
    Population city_popu;
    Population final_popu;
    feed_data(&CT);

    final_popu.pa = (Path*)malloc(POPULATION_SIZE * sizeof(Path));
    final_popu.no_of_path = POPULATION_SIZE;

    int path_per_slave = floor(POPULATION_SIZE / size);
    int cnt1 = path_per_slave;
    int cnt2 = path_per_slave;
    int generation_count1 = 0;
    int generation_count2 = 0;
    int generation_count3 = 0;
    if(rank == 0){
        create_Pop(&city_popu,path_per_slave,rank);
        cal_fit_Pop(&city_popu, &CT);
        // Sort all paths in current population, the path has biggest fitness, mean smallest length will be the first
        sort_fitness(&city_popu);
        // Calculate and update path probability:
        update_probability(&city_popu);
        printf("\n");
        //print_pop(&city_popu);
        while(1){
            generation_count1 ++;            
            order_crossover(&city_popu);
            mutation(&city_popu);
            cal_fit_Pop(&city_popu, &CT);
            sort_fitness(&city_popu);
            update_probability(&city_popu);

            // After a period, send the best individual to neighbor processor
            // So master will send to processor 1 and processor (size -1), and also receive from these processors
            if(generation_count1 % MIGRATION_PERIOD == 0  && generation_count1 > 0){
                //print_pop(&city_popu);
                MPI_Send(&city_popu.pa[0], sizeof(Path), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                MPI_Send(&city_popu.pa[0], sizeof(Path), MPI_CHAR, size-1, 0, MPI_COMM_WORLD);
                MPI_Recv(&city_popu.pa[path_per_slave-1], sizeof(Path), MPI_CHAR, 1, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                MPI_Recv(&city_popu.pa[path_per_slave-2], sizeof(Path), MPI_CHAR, size-1, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
            }
            if(generation_count1 == NUM_GENERATION){
                cnt2 = path_per_slave;
                for (i = 0; i < path_per_slave; i++) {
                    final_popu.pa[i] = city_popu.pa[i];
                }
                for(i=1; i<size; i++){
                    for(j=0; j<path_per_slave; j++){
                    MPI_Recv(&final_popu.pa[cnt2], sizeof(Path), MPI_CHAR, i, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                    cnt2++;
                    }
                }
                sort_fitness(&final_popu);
                print_pop(&final_popu);
                break;
            }

        }
    }
    else if (rank == size-1){
        create_Pop(&city_popu,path_per_slave,rank);
        cal_fit_Pop(&city_popu, &CT);
        // Sort all paths in current population, the path has biggest fitness, mean smallest length will be the first
        sort_fitness(&city_popu);
        // Calculate and update path probability:
        update_probability(&city_popu);
        printf("\n");
        while(1){
            generation_count2 ++;            
            order_crossover(&city_popu);
            mutation(&city_popu);
            cal_fit_Pop(&city_popu, &CT);
            sort_fitness(&city_popu);
            update_probability(&city_popu);
            // After a period, send the best individual to neighbor processor
            // So master will send to processor 1 and processor (size -1), and also receive from these processors
            if(generation_count2 % MIGRATION_PERIOD == 0  && generation_count2 > 0){
                MPI_Send(&city_popu.pa[0], sizeof(Path), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&city_popu.pa[0], sizeof(Path), MPI_CHAR, size -2, 0, MPI_COMM_WORLD);
                MPI_Recv(&city_popu.pa[path_per_slave-1], sizeof(Path), MPI_CHAR, 0, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                MPI_Recv(&city_popu.pa[path_per_slave-2], sizeof(Path), MPI_CHAR, size-2, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
            }
            if(generation_count2 == NUM_GENERATION){
                for(j=0; j<path_per_slave; j++){
                    MPI_Send(&city_popu.pa[j], sizeof(Path), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                }
                break;
           }
        }   
    }
    else {
        create_Pop(&city_popu,path_per_slave,rank);
        cal_fit_Pop(&city_popu, &CT);
        // Sort all paths in current population, the path has biggest fitness, mean smallest length will be the first
        sort_fitness(&city_popu);
        // Calculate and update path probability:
        update_probability(&city_popu);
        printf("\n");
        while(1){
            generation_count3 ++;            
            order_crossover(&city_popu);
            mutation(&city_popu);
            cal_fit_Pop(&city_popu, &CT);
            sort_fitness(&city_popu);
            update_probability(&city_popu);
            // After a period, send the best individual to neighbor processor
            // So master will send to processor 1 and processor (size -1), and also receive from these processors
            if(generation_count3 % MIGRATION_PERIOD == 0  && generation_count3 > 0){
                MPI_Send(&city_popu.pa[0], sizeof(Path), MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send(&city_popu.pa[0], sizeof(Path), MPI_CHAR, rank -1, 0, MPI_COMM_WORLD);
                MPI_Recv(&city_popu.pa[path_per_slave-1], sizeof(Path), MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                MPI_Recv(&city_popu.pa[path_per_slave-2], sizeof(Path), MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
            }
            if(generation_count3 == NUM_GENERATION){
                for(j=0; j<path_per_slave; j++){
                    MPI_Send(&city_popu.pa[j], sizeof(Path), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                }
                break;
            }
        }   
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    end_time = MPI_Wtime();
    if(rank == 0) { 
        printf("* Total parallel time = %f\n", end_time-start_time);
    }
    //After malloc, free the memory
    free(city_popu.pa);
    free(final_popu.pa);
    free(CT.city_p);    
    MPI_Finalize();
    return 0;
}
