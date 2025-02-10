#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#define N_FIELDS 2
#define T_FINAL 10.0 // Tempo final
#define dt 0.0001 // Discretização temporal

float 
    pmin = 0.001, // Densidade mínima espécie P
    pmax = 1.0, // Densidade máxima espécie P
    rh = 0.5, // Taxa de crescimento espécie H
    rp = 0.5, // Taxa de crescimento espécie P
    w1 = 1, // Peso competição espécie H
    w2 = 1, // Peso competição espécie P
    Dh = 0.7, // Coeficiente de difusão espécie H
    Dp = 0.3, // Coeficiente de difusão espécie P
    domain_size = 10.0; // Tamanho do domínio 

float source_points[1] = {-1.0};

enum CellType {
    H, 
    P
};

typedef struct TAxis {
    float min;
    float max; 
    float dim;
    int size;
    float delta;
    float *values;
} Axis; 

typedef struct TField {
    enum CellType cell_type;
    float *cur_values;
    float *next_values;
    Axis x_axis;
    Axis y_axis;
    char name[100];
    float *source_points;
    float *sink_points;
} Field;

float new_precision(float v, int i){
    return floorf(powf(10,i)*v)/powf(10,i);
}

int get_index(Field *field, int i, int j){
    return i + j*(field->x_axis.size);
}

float get_current_value(Field *field, int i, int j){
    return field->cur_values[get_index(field, i, j)];
}

float get_next_value(Field *field, int i, int j){
    return field->next_values[get_index(field, i, j)];
}

void set_current_value(Field *field, int i, int j, float v){
    field->cur_values[get_index(field, i, j)] = v;
}

void set_next_value(Field *field, int i, int j, float v){
    field->next_values[get_index(field, i, j)] = v;
}

inline float g(float w){
    return (1.0 - w) < 0.0? 0.0: (1.0 - w);
}

inline float get_sum(Field** fields, int x, int y){
    float sum = 0.0;
    sum = w1*get_current_value(fields[0], x, y) + w2*get_current_value(fields[1], x, y);
    /*for (int i = 0; i < N_FIELDS; i++){
        sum += get_current_value(fields[i],x, y);
    }*/
    return sum;
}

float chemotaxis(Field *p_field, Field *ch_field, int x, int y){
    float resx = 0, resy = 0, flux_left = 0, flux_right = 0;
    float p = get_current_value(p_field, x, y);
    float ch = get_current_value(ch_field, x, y);
    float dx = p_field->x_axis.delta, dy = p_field->y_axis.delta;

    if (x > 0){
        float pm1 = get_current_value(p_field, x - 1, y);
        float chm1 = get_current_value(ch_field, x - 1, y);
        if(ch -  chm1 > 0)        
            flux_left = -((ch - chm1) * pm1)/dx;
        else        
            flux_left = -((ch - chm1) * p)/dx;
    }
    if(x < p_field->x_axis.size - 1){
        float pp1 = get_current_value(p_field, x + 1, y);
        float chp1 = get_current_value(ch_field, x + 1, y);
        if(chp1 - ch > 0)        
            flux_right = (chp1 - ch)*p /dx;         
        else        
            flux_right = (chp1 - ch)*pp1 /dx;        
    }
    resx = (flux_left + flux_right)/dx;

    flux_left = 0, flux_right = 0;
    if (y > 0) {
        
        float pm1 = get_current_value(p_field, x, y - 1);
        float chm1 = get_current_value(ch_field, x, y - 1);

        if(ch - chm1 > 0)        
            flux_left = -(ch - chm1)*pm1 /dy;         
        else
            flux_left = -(ch - chm1)*p /dy;
    }
    if(y < p_field->y_axis.size - 1) {
        
        float pp1 = get_current_value(p_field, x, y + 1);
        float chp1 = get_current_value(ch_field, x, y + 1);

        if(chp1 - ch > 0)        
            flux_right = (chp1 - ch)*p /dy; 
        else 
            flux_right = (chp1 - ch)*pp1 /dy; 
    }
    resy = (flux_left + flux_right)/dy;

    return resx + resy; 
}

float diffusion(Field *field, int x, int y){
    float ux = get_current_value(field, x, y), uxp1, uxm1, diff_x, diff_y;
    float dx = field->x_axis.delta;
    float dy = field->y_axis.delta;

    if (x == 0) {
        uxp1 = get_current_value(field, x + 1, y);
        diff_x = (uxp1 - ux) /(dx*dx);
    }
    else if (x == (field->x_axis.size - 1)) {        
        uxm1 = get_current_value(field, x - 1, y);
        diff_x = (uxm1 - ux)/(dx*dx);
    }
    else {
        uxp1 = get_current_value(field, x + 1, y);   
        uxm1 = get_current_value(field, x - 1, y);
        diff_x = (uxp1 - 2*ux + uxm1)/(dx*dx);
    }
    
    if (field->y_axis.size == 1) {
        diff_y = 0;
    }
    else {
        float uy = get_current_value(field, x, y), uyp1, uym1;
        uyp1 = get_current_value(field, x, y + 1);
        uym1 = get_current_value(field, x, y - 1);
        if (y == 0) {
            diff_y = (uyp1 - uy) /(dy*dy);
        }
        else if (y == (field->y_axis.size - 1)) {
            diff_y = (uym1 - uy)/(dy*dy);
        }
        else {    
            diff_y = (uyp1 - 2*uy + uym1)/(dy*dy);
        }
    }
    
    return diff_x + diff_y;
}

float advection(float v, Field* field, int i, int j){
    float coeff = v*(1/field->x_axis.delta);
    float advec = 0;
    
    if (v > 0) {
        if (i > 0)
            advec = coeff *(get_current_value(field, i, j) - get_current_value(field, i - 1, j));
    }
    else if (v < 0) {
        if (i < field->x_axis.size - 1)
            advec = coeff *(get_current_value(field, i + 1, j) - get_current_value(field, i, j));
    }
    return advec;
}

void save_current_field(Field *field, float t) {
    char time[15];
    sprintf(time, "_%.3f.csv", t);
    char f_name[100];
    strcpy(f_name,field->name);
    strcat(f_name, time);
    FILE *fp = fopen(f_name, "w");
    fprintf(fp, "x,y,value\n");
    for (int j = 0; j < field->y_axis.size; j++) {
        for (int i = 0; i < field->x_axis.size; i++)
            fprintf(fp, "%.5f,%.5f,%.5f\n", field->x_axis.values[i], field->y_axis.values[j], get_current_value(field, i, j));
    }   
    fprintf(fp,"\n"); 
    fclose(fp);
}

void save_field(Field *field, float t) {
    char time[15];
    sprintf(time, "_%.3f.csv", t);
    char f_name[100];
    strcpy(f_name,field->name);
    strcat(f_name, time);
    FILE *fp = fopen(f_name, "w");
    fprintf(fp, "x,y,value\n");
    for (int j = 0; j < field->y_axis.size; j++) {            
        for (int i = 0; i < field->x_axis.size; i++)
            fprintf(fp, "%.5f,%.5f,%.5f\n", field->x_axis.values[i], field->y_axis.values[j], get_current_value(field, i, j));
    }   
    fprintf(fp,"\n"); 
    fclose(fp);
}

Axis create_axis(float min, float max, float dx){
    Axis axis;
    axis.min = min;
    axis.max = max;
    axis.dim = (max - min);
    axis.delta = dx;
    axis.size = (axis.dim/dx); 
    printf("size: %d\n", axis.size);

    axis.values = (float*) calloc(axis.size, sizeof(float));
    
    int i = 0;
    for (float x = 0; x < axis.dim; x += axis.delta) {
        axis.values[i] = new_precision(x, 2.0);
        i++;
    }
    return axis;
}

Field* create_field(enum CellType cell_type, char *name, Axis x_axis, Axis y_axis){
    Field *field = (Field*) calloc(1, sizeof(Field));
    strcpy(field->name, name);
    field->cell_type = cell_type;
    field->x_axis = x_axis;
    field->y_axis = y_axis;

    int field_size = x_axis.size*y_axis.size;
    field->cur_values = (float*) calloc(field_size, sizeof(float));
    field->next_values = (float*) calloc(field_size, sizeof(float));
    return field;
}

Field** create_fields(enum CellType cell_type[], char *name[], Axis x_axis, Axis y_axis){    
    Field** fields = (Field**) calloc(N_FIELDS, sizeof(Field*));
    for (int i = 0; i < N_FIELDS; i++){
        fields[i] = create_field(cell_type[i], name[i], x_axis, y_axis);
    }
    return fields;
}

int binary_search(float x, float *array, int low, int high) {    
    while (low <= high) {
        int mid = (low + high) / 2;

        if (fabsf(array[mid] - x) < powf(10.0,-3.0)) {
            return mid;
        }

        if (x > array[mid])
            low = mid + 1;

        else
            high = mid - 1;
    }
    return -1;
}

inline float derivation_x(Field* field, int i, int j){
    float central_diff = 0;
    if (i == 0){
        central_diff = (get_current_value(field, i + 1, j) - get_current_value(field, i, j))/ (field->x_axis.delta);
    }
    else if (i == field->x_axis.size - 1){
        central_diff = (get_current_value(field, i, j) - get_current_value(field, i - 1, j))/ (field->x_axis.delta);
    }
    else {
        central_diff = (get_current_value(field, i + 1, j) - get_current_value(field, i - 1, j))/ (2*field->x_axis.delta);
    }
    return central_diff;
}

float get_fields_sum_at_pos(Field** fields, int i, int j, int indexes[], int len){
    float sum  = 0.0;
    for (int i = 0; i < len; i++){
        sum += get_current_value(fields[indexes[i]], i, j);
    }
    return sum;
}


inline int get_index_coord(float x, float *array, int len){ 
    return binary_search(x, array, 0, len - 1); 
}

float Pmin(Field* field, float x[], int len, int i, int j){
    if (x[0] == -1) 
        return pmin;
    
    for (int k = 0; k < len; k++){
        if ((x[k]/field->x_axis.delta) == i)
            return pmin;   
    }
    return 0.0;
}

float Pmax(Field* field, float x[], int len, int i, int j){
    if (x[0] == -1) 
        return pmax;
    
    for (int k = 0; k < len; k++){
        if ((x[k]/field->x_axis.delta) == i)
            return pmax;   
    }

    return 0.0;
}

inline float g_right_border_avg(Field** fields, int x, int y){
    return g((get_sum(fields,x+1,y)+get_sum(fields,x,y))/2.0 );
}

inline float g_up_border_avg(Field** fields, int x, int y){
    return g((get_sum(fields,x,y+1)+get_sum(fields,x,y))/2.0 );
}

inline float g_left_border_avg(Field** fields, int x, int y){
    return g((get_sum(fields,x,y)+get_sum(fields,x-1,y))/2.0 );
}

inline float g_down_border_avg(Field** fields, int x, int y){
    return g((get_sum(fields,x,y)+get_sum(fields,x,y-1))/2.0 );
}

inline float u_right_border(Field* u, int x, int y){
    return get_current_value(u, x+1, y) - get_current_value(u, x, y);
}

inline float u_up_border(Field* u, int x, int y){
    return get_current_value(u, x, y+1) - get_current_value(u, x, y);
}

inline float u_left_border(Field* u, int x, int y){
    return get_current_value(u, x, y) - get_current_value(u, x-1, y);
}

inline float u_down_border(Field* u, int x, int y){
    return get_current_value(u, x, y) - get_current_value(u, x, y-1);
}

inline float local_avg_x(Field** fields, Field* u, int x, int y){
    float dx = u->x_axis.delta;
    if(x == 0)
        return (g_right_border_avg(fields, x, y)*u_right_border(u, x, y)) /(dx*dx);
    else if(x == u->x_axis.size-1)
        return - (g_left_border_avg(fields, x, y)*u_left_border(u, x, y)) /(dx*dx);
    else
        return (g_right_border_avg(fields, x, y)*u_right_border(u, x, y) - g_left_border_avg(fields, x, y)*u_left_border(u, x, y))/(dx*dx);
}

inline float local_avg_y(Field** fields, Field* u, int x, int y){
    float dy = u->y_axis.delta;
    if(y == 0)
        return (g_up_border_avg(fields, x, y)*u_up_border(u, x, y)) /(dy*dy);
    else if(y == u->y_axis.size - 1)
        return - (g_down_border_avg(fields, x, y)*(u_down_border(u, x, y))) /(dy*dy);
    else
        return (g_up_border_avg(fields, x, y)*u_up_border(u, x, y) -  g_down_border_avg(fields, x, y)*u_down_border(u, x, y) ) /(dy*dy);
}

inline float local_average(Field** fields, Field* u, int x, int y){
    return local_avg_x(fields, u, x, y) + local_avg_y(fields, u, x, y);
}

int main(){    
    int nsteps = (int)((T_FINAL+dt)/dt);
    printf("steps: %d\n", nsteps);   
    int interval = (int)(T_FINAL/dt)/10;
    
    FILE *tFile = fopen("results/t.csv", "w");
    for (float t = 0; t < T_FINAL+dt; t += dt)
        fprintf(tFile, "%.6lf\n", t);
    fclose(tFile);    

    enum CellType field_types[N_FIELDS] = {H, P};
    char *field_names[N_FIELDS] = {"results/H", "results/P"};

    float axis_min = 0.0, axis_max = domain_size, axis_delta = 0.1;
    Axis x_axis = create_axis(axis_min, axis_max, axis_delta);
    Axis y_axis = create_axis(axis_min, axis_max, axis_delta);

    Field** fields = create_fields(field_types, field_names, x_axis, y_axis);
    Field* h_field = fields[0];
    Field* p_field = fields[1];
    
    FILE *xFile = fopen("results/x.csv", "w"); 
    for (float x = h_field->x_axis.min; x <= h_field->x_axis.max ; x +=  h_field->x_axis.delta) {
        fprintf(xFile, "%.2lf\n", x);        
    }
    fclose(xFile); 

    FILE *yFile = fopen("results/y.csv", "w"); 
    for (float y = h_field->y_axis.min; y <= h_field->y_axis.max ; y +=  h_field->y_axis.delta) {
        fprintf(yFile, "%.2lf\n", y);
    }
    fclose(yFile); 

    int x_size = h_field->x_axis.size;
    int y_size = h_field->y_axis.size;
    int num = 100, k = 0;
    srand(time(NULL));
    //condicao inicial 
    while (k < num){
        int pos_x = rand() % x_size;
        int pos_y = rand() % y_size;
        set_current_value(h_field, pos_x, pos_y, 0.2);
        k++;
    }

    k = 0;
    while (k < num){
        int pos_x = rand() % x_size;
        int pos_y = rand() % y_size;
        set_current_value(p_field, pos_x, pos_y, 0.2);
        k++;
    }
   
    for (int i =0; i < N_FIELDS; i++)
        save_current_field(fields[i], 0);    

    float t = 0;
    for(int step = 1; step <= nsteps; ++step) {
                    
        for (int j = 0; j < y_size; j++) {

            #pragma omp parallel for num_threads(6) 
            for (int i = 0; i < x_size; i++) {

                float h = get_current_value(h_field, i, j);
                float p = get_current_value(p_field, i, j);
                //testa regiao 
                if (i < x_size/2) {
                    w1 = 1;
                    w2 = 1;
                }
                else {
                    w1 = 1;
                    w2 = 4;
                }
                
                float h_next = rh*h*g(get_sum(fields, i, j)) + Dh*local_average(fields, h_field, i, j);

                float p_next = rp*p*g(get_sum(fields, i, j)) + Dp*local_average(fields, p_field, i, j);

                set_next_value(h_field, i, j,  h + (h_next * dt));                
                set_next_value(p_field, i, j, p + (p_next * dt));
                
            }
        }

        if (step % interval == 0) {
            for (int i =0; i < N_FIELDS; i++)
                save_field(fields[i], t);      
            printf("saving time %f\n", t);
        }

        for (int i =0; i < N_FIELDS; i++) {
            float *aux = fields[i]->cur_values;
            fields[i]->cur_values = fields[i]->next_values;
            fields[i]->next_values = aux;
        }

        t += dt;  
    }

    for (int i = 0; i < N_FIELDS; i++) {
        free(fields[i]->cur_values);
        free(fields[i]->next_values);
        free(fields[i]);
    }
    free(fields);
    return 0; 
}
