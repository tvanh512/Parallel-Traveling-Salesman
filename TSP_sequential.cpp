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
    float total_fit;
};


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

int create_Path(Path *path_p){
    // First, create a temporary list of cities
    int j;
    for(j=0; j<NUM_CITIES; j++){
        path_p->citylist[j] = j; 
    }
    // Then randomly select a position, then swap it with the current position.
    int i;    
    for (i = NUM_CITIES-1; i > 0; i--) {
        int j = (unsigned int) (drand48()*(i+1));
        int t = path_p->citylist[j];
        path_p->citylist[j] = path_p->citylist[i];
        path_p->citylist[i] = t;
    }
}

int create_Pop(Population *popu){
    popu->pa = (Path*)malloc(POPULATION_SIZE * sizeof(Path));
    int i;
    for(i=0; i<POPULATION_SIZE; i++){
        create_Path(&popu->pa[i]);
    }
}

int cal_fit_Pop(Population *popul, Country *count){
    int i;
    float total_Fit = 0.0;

    for(i=0; i<POPULATION_SIZE; i++){
        popul->pa[i].distance = cal_Fit(&popul->pa[i].citylist[0], count);
        popul->pa[i].fitness = 1000000/popul->pa[i].distance;
        total_Fit = total_Fit + popul->pa[i].fitness;
    }
    popul->total_fit = total_Fit;
}

int cmpfunc (const void *a, const void *b){
  Path * P1 = (Path *)a;
  Path * P2 = (Path *)b;
  // The return type is int so multiple by a big number help when the distance differ by a small fraction after dot
  return (P1->distance - P2->distance);
}

float sort_fitness(Population *popu){
    qsort(popu->pa, POPULATION_SIZE, sizeof(Path), cmpfunc);
}

void update_probability(Population *popu){
    int i;
    for (i=0; i < POPULATION_SIZE; i++){
        popu->pa[i].probability = popu->pa[i].fitness/popu->total_fit;
    }
}

void print_pop(Population *popu){
    printf("The current best path distance is: \n");
    printf("%f \n",popu->pa[0].distance);   
    printf("-------------------- \n");
}

float Roulette_Wheel_Selection(Population *popu){
    float sum;
    sum = 0;
    float r = drand48();
    int i;
    for (i=0; i < POPULATION_SIZE; i++){
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
    int mother;int father;
    
    mother = Roulette_Wheel_Selection(popu);
    father = Roulette_Wheel_Selection(popu);
    while (father == mother){
        father = Roulette_Wheel_Selection(popu);
    }
    int OX ;int temp;int i;    
    
    // Generate random number OX in range 1 to NUM_CITIES - 1, run the random several times.
    for (i=0; i < 3; i++){
        OX = (NUM_CITIES -1)*(double)rand()/RAND_MAX + 1;
    }    
    int son[NUM_CITIES] = {0};
    int daughter[NUM_CITIES] = {0};
    inherite_child(OX,son,&popu->pa[father],&popu->pa[mother]);
    inherite_child(OX,daughter,&popu->pa[mother],&popu->pa[father]);
    for (i=0; i<NUM_CITIES; i++){
            popu->pa[POPULATION_SIZE-1].citylist[i] = son[i];
            popu->pa[POPULATION_SIZE-2].citylist[i] = daughter[i];
        }
}

float mutation(Population *popu){
    // rand()%100 create a random number between 0 to 99
    if(rand()%100 < MUTATION_RATIO){
        int r;
        // r is a random number between 2 and POPULATION_SIZE - 1, because we don't want to mutate the 2 best paths
        r = 2+rand() % (POPULATION_SIZE-2) ;
        int i; 
        for (i = NUM_CITIES-1; i > 0; i--) {
            int j = (unsigned int) (drand48()*(i+1));
            int t = popu->pa[r].citylist[j];
            popu->pa[r].citylist[j] = popu->pa[r].citylist[i];
            popu->pa[r].citylist[i] = t;
        }
    }
}
int main()
{
    double start_time = time(0);
    Country CT;
    Population city_popu;
    feed_data(&CT);
    // Initialization
    create_Pop(&city_popu);
    int i;float dist ;float fit ;       

    for (i = 0; i < NUM_GENERATION; i++) {
        // Calculate the fitness of every path, and the total fitness of the population
        cal_fit_Pop(&city_popu, &CT);
        // Sort all paths in current population, the path has biggest fitness, mean smallest length will be the first
        sort_fitness(&city_popu);
        // Calculate and update path probability:
        update_probability(&city_popu);
        /*if (i%500==0){
            print_pop(&city_popu);
        }*/
        order_crossover(&city_popu);
        mutation(&city_popu);
    }

    cal_fit_Pop(&city_popu, &CT);
    sort_fitness(&city_popu);
    update_probability(&city_popu);
    print_pop(&city_popu);
    double end_time = time(0);;
    printf("* Total sequential time = %f\n", end_time-start_time);
    return 0;
}
