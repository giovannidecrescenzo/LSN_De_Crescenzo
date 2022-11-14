#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm> 
#include <numeric>
#include "population.h"
#include "../Random/random.h"

using namespace std;

Population :: Population(int cities, int num_chromo, Random rnd, int t){
    n_cities = cities;
    n_chromo = num_chromo;
    class_rnd = rnd;
    type = t;
    Init_cities();
    Init_pop();
    Fill_fit();
}

Population :: Population(int cities, int num_chromo, Random rnd, int t, string filename){
    n_cities = cities;
    n_chromo = num_chromo;
    class_rnd = rnd;
    type = t;
    Init_cities_input(filename);
    Init_pop();
    Fill_fit();
}

Population :: ~Population(){}

double Population :: Fit(vector<int> chromo){
    if(type == 0.){
        double sum = 0.;
        for(int i = 1; i < chromo.size(); i++){
            sum += abs(2*sin((cities[chromo[i]]-cities[chromo[i-1]])/2));
        }
        sum += abs(2*sin(-(cities[chromo.size()-1]-cities[0])/2));
        return sum;
    }
    else{
        double sum = 0.;
        for(int i = 1; i < chromo.size(); i++){
            sum += sqrt((pow(sq_cities[chromo[i]][0]-sq_cities[chromo[i-1]][0],2) + pow(sq_cities[chromo[i]][1]-sq_cities[chromo[i-1]][1],2)));

        }

        sum += sqrt((pow(sq_cities[chromo[chromo.size()-1]][0]-sq_cities[chromo[0]][0],2)+pow(sq_cities[chromo[chromo.size()-1]][1]-sq_cities[chromo[0]][1],2)));
        return (sum);        
    }
}

void Population :: Init_cities(void){
    cities.clear();
    sq_cities.clear();
    if(type == 0){
        cities.push_back(0.);
        for(int i = 0; i < n_cities-1; i++){
            cities.push_back(class_rnd.Rannyu()*2*M_PI);
        }
    }
    else {
        double a = class_rnd.Rannyu();
        double b = class_rnd.Rannyu();
        vector<double> v(2);
        v[0] = a;
        v[1] = b;
        sq_cities.push_back(v);
        for(int i = 0; i < n_cities-1; i++){
            double x = class_rnd.Rannyu();
            double y = class_rnd.Rannyu();
            vector<double> v1(2);
            v1[0] = x;
            v1[1] = y;
            sq_cities.push_back(v1);            
        }
    }
}

void Population :: Init_cities_input(string filename){
    ifstream in;
    in.open(filename);
    for(int i = 0; i < n_cities; i++){
        vector<double> v(2);
        v[0] = 0;
        v[1] = 0;
        sq_cities.push_back(v);
    } 
    int i = 0; 
    while(i < n_cities){
        in >> sq_cities[i][0];
        in >> sq_cities[i][1];
        i++;
    }
    in.close();
}

void Population :: Init_pop(void){
    population.clear();
    vector<int> first;
    for(int i = 0; i < n_cities; i++){
        first.push_back(i);
    }
    population.push_back(first);
    for(int j = 0; j < n_chromo-1; j++){
        vector<int> temp;
        temp = population.back();
        for(int i = 1; i < 10; i++){
            rand_swap(temp);
        }
        population.push_back(temp);
    }
}

void Population :: Fill_fit(void){
    fitness.clear();
    for(int i = 0; i < n_chromo; i++){
        fitness.push_back(Fit(population[i]));
    }
}

void Population :: mixing(void){
    int i = int(class_rnd.Rannyu(0,n_chromo)); //decide with coin flip the chromosome to be modified
    double probabilities[6] = {1./7, 2./7, 3./7, 4./7, 4.3/7, 6./7};
    double x_1 = class_rnd.Rannyu();
    if(class_rnd.Rannyu() < 0.2){
        if(x_1 < probabilities[0]){
            rand_swap(population[i]);
            fitness[i] = Fit(population[i]);
        }
        if(x_1 > probabilities[0] and x_1 < probabilities[1]){
            translation(population[i]);
            fitness[i] = Fit(population[i]);
        }
        if(x_1 > probabilities[1] and x_1 < probabilities[2]){
            anti_translation(population[i]);
            fitness[i] = Fit(population[i]);
        }
        if(x_1 > probabilities[2] and x_1 < probabilities[3]){
            multi_swap(population[i]);
            fitness[i] = Fit(population[i]);
        }
        if(x_1 > probabilities[3] and x_1 < probabilities[4]){
            cont_Swap(population[i]);
            fitness[i] = Fit(population[i]);
        }
        if(x_1 > probabilities[4] and x_1 < probabilities[5]){
            inversion_swap(population[i]);
            fitness[i] = Fit(population[i]);
        }
        if(x_1 > probabilities[5]){
            cont_swap(population[i]);
            fitness[i] = Fit(population[i]);
        }
    }
    sort_population();
}




void Population :: Mutation(void){
    sort_population();
    if(class_rnd.Rannyu() < 0.95){
        //PARAMETERS TO BE REGULATED:
        //probabilities ratios of mutation over children
        //{0.5/7, 3/7, 5/7, 5.3/7, 5.6/7, 6.8/7};
        double probabilities[6] = {1./7, 2./7, 3./7, 3.5/7, 4.3/7, 6./7};
        //first interval 0,prob[0] regulates call of random swap of cities
        //second interval regulates translation of 1 position
        //third interval regulates anti_translation of 1 position
        //fourth interval regulates switch of two chunks of chromosome
        //fifth interval reglates  switch of two contigous cities of chromosome
        //sixth interval regulates mirroring
        //last interval prob[5] - 1 regulates switch of two contigous chunks of chromosome
        // p to regulate the weight of fittest parents
        double p = 4.5;

        int j_1, j_2;
        j_1 = int(n_chromo*pow(class_rnd.Rannyu(),p)); 
        do{
            j_2 = int(n_chromo*pow(class_rnd.Rannyu(),p)); 
        } while(j_1 == j_2);
        sort_population();
        vector<int> parent1, parent2;
        parent1 = population[j_1];
        parent2 = population[j_2];
        while (parent1==parent2){
            j_2 = int(n_chromo*pow(class_rnd.Rannyu(),p));
            parent2 = population[j_2];
        }

        int cut = int(class_rnd.Rannyu()*(n_cities-2))+1;
        int num = int(class_rnd.Rannyu()*(n_cities-cut));

         
        vector<vector<int> > children = crossover(parent1, parent2, cut, num);
        double x_1 = class_rnd.Rannyu();
        double x_2 = class_rnd.Rannyu();
        bool prob = false;
        if(class_rnd.Rannyu() < 0.1) prob = true; //probability for mutation on children

        if(x_1 < probabilities[0] and prob){
            rand_swap(children[0]);
        }
        if(x_2 < probabilities[0] and prob){
            rand_swap(children[1]);
        }
        if(x_1 > probabilities[0] and x_1 < probabilities[1] and prob){
            translation(children[0]);
        }
        if(x_2 > probabilities[0] and x_2 < probabilities[1] and prob){
            translation(children[1]);
        }

        if(x_1 > probabilities[1] and x_1 < probabilities[2] and prob){
            anti_translation(children[0]);
        }

        if(x_2 > probabilities[1] and x_2 < probabilities[2] and prob){
            anti_translation(children[1]);
        }

        if(x_1 > probabilities[2] and x_1 < probabilities[3] and prob){
            multi_swap(children[0]);
        }

        if(x_2 > probabilities[2] and x_2 < probabilities[3] and prob){
            multi_swap(children[1]);
        }

        if(x_1 > probabilities[3] and x_1 < probabilities[4] and prob){
            cont_Swap(children[0]);
        }

        if(x_2 > probabilities[3] and x_2 < probabilities[4] and prob){
            cont_Swap(children[1]);
        }
      
        if(x_1 > probabilities[4] and x_1 < probabilities[5] and prob){
            inversion_swap(children[0]);
        }

        if(x_2 > probabilities[4] and x_2 < probabilities[5] and prob){
            inversion_swap(children[1]);
        }

        if(x_1 > probabilities[5] and prob){
            cont_swap(children[0]);
        }

        if(x_2 > probabilities[5] and prob){
            cont_swap(children[1]);
        }
        sort_population();
        double fitness_1 = Fit(children[0]);
        double fitness_2 = Fit(children[1]);
        bool first = true;
        bool second = true;
        for(int j = 0; j < n_chromo; j++){
            if(children[0] == population[j]){
                first = false;
            }
            if(children[1] == population[j]){
                second = false;
            }
        }
        if(fitness_1 < fitness_2 and fitness_2 <= fitness.back() and first == true and second == true){
            population[n_chromo-2] = children[0];
            population[n_chromo-1] = children[1];
            fitness[n_chromo-2] = fitness_1;
            fitness[n_chromo-1] = fitness_2;            
        }
        if(fitness_2 < fitness_1 and fitness_1 <= fitness.back() and first == true and second == true){
            population[n_chromo-2] = children[0];
            population[n_chromo-1] = children[1];
            fitness[n_chromo-2] = fitness_1;
            fitness[n_chromo-1] = fitness_2;            
        }
        if(fitness_1 < fitness_2 and fitness_2 > fitness.back() and fitness_1 < fitness.back() and first == true and second == false){
            population[n_chromo-1] = children[0];
            fitness[n_chromo-1] = fitness_1;
        }
        if(fitness_2 < fitness_1 and fitness_1 > fitness.back() and fitness_2 < fitness.back() and second == true and first == false){
            population[n_chromo-1] = children[1];
            fitness[n_chromo-1] = fitness_2;            
        }
        sort_population();
    }
}


vector<size_t> Population :: sort_indices(vector<double> v) {


    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

void Population :: sort_population(void){
    vector<size_t> idx = sort_indices(fitness);
    vector<vector<int> > pop1 = population;
    vector <double> fit1 = fitness;
    for(int i = 0; i < fitness.size(); i++){
        population[i] = pop1[idx[i]];
        fitness[i] = fit1[idx[i]];
    }
}

vector<vector<int> > Population :: crossover(vector<int> a, vector<int> b, int cut, int num){
    vector<int> temp = a;
    vector<int> copy_a = a;
    vector<int> copy_b = b;
    vector<int> indices_a;
    vector<int> indices_b;
    vector<vector<int> > output;
    for(int i = cut; i < cut+num; i++){
        indices_a.push_back(getIndex(copy_b , copy_a[i]));
            
    }
    sort(indices_a.begin(), indices_a.end());
    for(int i = cut; i < cut+num; i++){
        copy_a[i] = copy_b[indices_a.front()];
        indices_a.erase(indices_a.begin());
    }
    output.push_back(copy_a);

    for(int i = cut; i < cut+num; i++){
        indices_b.push_back(getIndex(temp,b[i]));
    }
    sort(indices_b.begin(), indices_b.end());
    for(int i = cut; i < cut+num; i++){
        copy_b[i] = temp[indices_b.front()];
        indices_b.erase(indices_b.begin());
    }
    output.push_back(copy_b);
   
    return output;
}

void Population :: check(int pos){
    if(population[pos][0] != 0){
        cout << "ERROR AT POSITION " << pos << ": THERES NOT THE STARTING CITY IN FIRST POS" << endl;
    }
    for(int i = 0; i < n_cities; i++){
        for(int j = 0; j < n_cities and j != i; j++){
            if(population[pos][i] == population[pos][j]){
                cout << "ERROR IN POSITION " << pos << ": 2 CITIES ARE THE SAME" << endl;
            }
        }
    }
}

//POSSIBLE MUTATIONS

//rand_swap: swap order of two random cities
void Population :: rand_swap(vector<int> &chromo){
    int a = int(class_rnd.Rannyu()*(n_cities-1))+1;
    int b = int(class_rnd.Rannyu()*(n_cities-1))+1;
    vector<int> temp = chromo;
    chromo[a] = temp[b];
    chromo[b] = temp[a];
}
//cont_swap: swap two contigous cities
void Population :: cont_Swap(vector<int> &chromo) {
    int cut = int(class_rnd.Rannyu(1,n_cities-2));
    vector<int> temp = chromo;
    chromo[cut] = temp[cut+1];
    chromo[cut + 1] = temp[cut];
}
//multi_swap: takes #num cities starting at cut and switches them with #num cities starting from cut1
void Population :: multi_swap(vector<int> &chromo){
    int cut = int(class_rnd.Rannyu(1, n_cities-1));
    int num = int(class_rnd.Rannyu()*int((n_cities-cut)/2))+1;
    int cut1;
    if(cut + num == n_cities - num - 1) cut1 = cut + num;

    else cut1 = int(class_rnd.Rannyu(cut + num, n_cities - num + 1));
    vector<int> temp = chromo;
    for(int i = 0; i < num; i++){
        chromo[i+cut] = temp[i+cut1];
        chromo[i+cut1] = temp[i+cut];
    }

}

void Population :: cont_swap(vector<int> &chromo){
    int cut = int(class_rnd.Rannyu()*(n_cities-2))+1;
    int num = int(class_rnd.Rannyu()*((n_cities-cut)/2));
    vector<int> temp = chromo;
    for(int i = cut; i < cut + num; i++){
        chromo[i] = temp[num+i];
        chromo[num+i] = temp[i];
    }

}


void Population :: translation(vector<int> &chromo){
    int cut = int(class_rnd.Rannyu()*(n_cities-2))+1;
    int num = int(class_rnd.Rannyu()*(n_cities-cut-1))+1;
    vector<int> temp = chromo;
    for(int i = cut; i < cut+num-1; i++){
        chromo[i] = temp[i+1];
    }
    chromo[cut+num-1] = temp[cut];
}

void Population :: anti_translation(vector<int> &chromo){
    int cut = int(class_rnd.Rannyu(1,n_cities-1));
    int num = int(class_rnd.Rannyu()*(n_cities-cut-1))+1;
    vector<int> temp = chromo;
    for(int i = cut+1; i <= cut+num; i++){
        chromo[i] = temp[i-1];
    }
    chromo[cut] = temp[cut+num];
}

//inversion_swap: mirror image of #num cities starting form cut
void Population :: inversion_swap(vector<int> &chromo){
    int cut = int(class_rnd.Rannyu()*(n_cities-1))+1;
    int num = int(class_rnd.Rannyu()*(n_cities-cut));
    vector<int> temp = chromo;
    for(int i = cut; i < int((cut + num)/2); i++){
        chromo[i] = temp[num+cut+1-i];
        chromo[num+cut+1-i] = temp[i];
    }
}

int Population :: getIndex(vector<int> v, int K)
{
    auto it = find(v.begin(), v.end(), K);
  
    // If element was found
    if (it != v.end()) 
    {
        // calculating the index
        // of K
        int index = it - v.begin();
        return index;
    }
    else {
        // If the element is not
        // present in the vector
        cout << "ELEMENT NOT PRESENT IN VECTOR" << endl;
        return -1;
    }
}
