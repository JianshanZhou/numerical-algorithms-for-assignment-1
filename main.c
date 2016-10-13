/* information
Copyright (C) Thu Oct 11 10:11:08 2016  Jianshan Zhou

Contact: zhoujianshan@buaa.edu.cn	jianshanzhou@foxmail.com

Website: <https://github.com/JianshanZhou>

This program is free software: you can redistribute
 it and/or modify it under the terms of
 the GNU General Public License as published
 by the Free Software Foundation,
 either version 3 of the License,
 or (at your option) any later version.

This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program.
 If not, see <http://www.gnu.org/licenses/>.

This is the main file where the main function and other sub-functions are provided
 to analyze the 1-st issue presented in the Numerical Analysis class.
*/


//import some necessary basic libs
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//define some macros
#define N 501
#define S 2
#define R 2


//tolerant accuracy
const double epsilon = 1.0e-12;
const double minU = 1.0, maxU = 3.0;


// declarations of all the sub-functions
int Initiate_Matrix(double *A);
int Initiate_CompressedMatrix(double **matrix, double *array);
double IndexMapping(double *arr, int i, int j);
int Exchange(double *Value1, double *Value2);
int Determine_Max_Min(double *Value1, double *Value2);
int min_two_values(int a, int b);
int max_three_values(int a, int b, int c);
int Doolittle_LU_composition(double **matrix, double *det_value);
double Sum_up(double **matrix, int t1, int t2, int k, int j);
double Sum_up2(double ** matrix, int t1, int t2, int k, int i);
double Norm(double *array);
int Multiply_matrix_and_vector(double *array, double *vector, double * output);
int Initiate_vector(double *array, double a, double b);
int Power_algorithm(double *array, double *eigVector, double *eigValue, double shift);
int inPower_algorithm(double *array, double *eigVector, double *eigValue, double shift);
double innerProduc(double *arr1, double *arr2);
int Multiply_matrix_and_vector2(double *array, double *vector, double *output, double shift);
double Sum_up3(double **matrix,double *b,int i,int t_lower,int t_upper);
double Sum_up4(double **matrix, double *solution, int i, int t_lower, int t_upper);
int Solve_linear_equation(double **matrix, double *b, double *solution);
int Initiate_CompressedMatrix_shift(double **matrix, double *array, double shift);
int write_data_to_txt(FILE *filepointer, double *array, int n);
int Analyze(void);
int problems(double * results);

//some testing functions
int test(void);

int main()
{
    //test();//test some functions
    printf("The whole program to solve the following three problems is developed by:\n");
    printf("Zhou Jianshan BY1613123\n");
    printf("Email: jianshanzhou@foxmail.com\n");
    printf("Web: https://github.com/JianshanZhou\n");
    printf("Please feel free to contact me if you have any questions!\n");

    /*Initialize the matrix A that is formulated as a 1-D array*/
    double array[N+2] = {0.0};
    Initiate_Matrix(array);

    /*Problem-1: solve for the maximum and the minimum eigenvalues of the matrix A*/
    double lambda1 = 0.0, lambda501 = 0.0, lambda_s = 0.0;
    double shift = 0.0;
    double eigVector[N] = {1.0};
    //solve for the first eig. value
    Power_algorithm(array, eigVector, &lambda1, shift);
    //printf("The first eig. value obtained is %.11e\n",lambda1);
    //let the shift to be the first lambda_1
    shift = lambda1;
    //solve for the second eig. value
    Power_algorithm(array, eigVector, &lambda501, shift);
    lambda501 = lambda501 + shift;
    //printf("The second eig. value obtained is %.11e\n",lambda501);
    //let lambda1 to be the minimum one while lambda501 the maximum
    if(lambda1>lambda501)
    {
        Exchange(&lambda1,&lambda501);
    }
    printf("The minimum eigenvalue is: lambda1 = %.11e\n",lambda1);
    printf("The maximum eigenvalue is: lambda501 = %.11e\n",lambda501);
    //solve for lambda_s by using the inverse Power Alg.
    shift = 0.0;
    inPower_algorithm(array, eigVector, &lambda_s, shift);
    printf("The min-abs eigenvalue is: lambda_s = %.11e\n",lambda_s);

    /*Problem-2: solve for lambda_i(k), k= 1,2,...,39*/
    int k_num = 39, k;
    double lambda_arr[39] = {1.0};
    for(k=1;k<=k_num;k++)
    {
        shift = lambda1 + k*((lambda501-lambda1)/40.0);
        inPower_algorithm(array, eigVector, &lambda_arr[k-1], shift);
        lambda_arr[k-1] = lambda_arr[k-1] + shift;
        printf("The eigenvalue lambda_i(%d) is %.11e closed to %.11e\n",k,lambda_arr[k-1],shift);
    }

    /*Problem-3: solve the 2-norm-based condition of the matrix A and its determinant det(A)*/
    double cond = fabs(lambda1/lambda_s);
    printf("The 2-norm-based condition of this matrix A is cond(A) = %.11e\n", cond);

    double det_value = 0.0;
    shift = 0.0;
    double matrix[R+S+1][N] = {0.0};
    Initiate_CompressedMatrix_shift(matrix, array, shift);
    //step-3: do Doolittle LU decomposition
    Doolittle_LU_composition(matrix, &det_value);//A=LU
    printf("The determinant of this matrix A is det(A) = %.11e\n", det_value);


    /*Analyze the algorithms by Monte Carlo simulations.*/
    //Analyze();

    return 0;
}



// define all of the sub-functions here
/*This function analyzes the algorithms developed here.*/
int Analyze(void)
{

//Analyze the impact of different randomly initialized u(0) on the numerical results
    int sim_num = 100;
    double tempS = 0.0, tempS2 = 0.0;
    double array_results[100][3+39+2] = {0.0};
    double avg_re[3+39+2] = {0.0};
    double std_re[3+39+2] = {0.0};
    int sim_flag = 0;
    for(sim_flag=1;sim_flag<=sim_num;sim_flag++)
    {
        printf("Running the %d-th Monte Carlo simulation!\n",sim_flag);
        problems(&array_results[sim_flag-1]);//the pointer to the first element's address of each row
    }
    //average all the results
    int k;
    for(k=0;k<(3+39+2);k++)
    {
        tempS = 0.0;
        for(sim_flag=0;sim_flag<sim_num;sim_flag++)
        {
            tempS += array_results[sim_flag][k];
        }
        avg_re[k] = tempS/(sim_num*1.0);
        tempS2 = 0.0;
        for(sim_flag=0;sim_flag<sim_num;sim_flag++)
        {
            if(fabs((array_results[sim_flag][k]-avg_re[k])*(array_results[sim_flag][k]-avg_re[k]))<=epsilon)
            {
                tempS +=epsilon;
            }
            else
            {
                tempS += (array_results[sim_flag][k]-avg_re[k])*(array_results[sim_flag][k]-avg_re[k]);
            }
        }
        if(sqrt(tempS/(sim_num))<=epsilon)
        {
            std_re[k] = epsilon;
        }
        else
        {
            std_re[k] = sqrt(tempS/(sim_num));
        }
    }

    //export the results into the external file *.txt
    /*get the file path*/
    char data_file_name1[] = "C:\\Users\\zhoujianshan\\OneDrive\\ÎÄµµ\\C_projects\\zhoujianshan_by1613123_project1\\avg_results.txt";
    char data_file_name2[] = "C:\\Users\\zhoujianshan\\OneDrive\\ÎÄµµ\\C_projects\\zhoujianshan_by1613123_project1\\std_results.txt";

    FILE *pWrite;
    if((pWrite=fopen(data_file_name1,"w"))==NULL){
        printf("Errors occur in creating a txt file pointer!\n");
        exit(1);
    }
    /*write the avg data to a txt file*/
    write_data_to_txt(pWrite, avg_re, 3+39+2);

    if((pWrite=fopen(data_file_name2,"w"))==NULL){
        printf("Errors occur in creating a txt file pointer!\n");
        exit(1);
    }
    /*write the std data to a txt file*/
    write_data_to_txt(pWrite, std_re, 3+39+2);

    return 0;
}

/*solve the problem-1-2-3*/
int problems(double * results)
{

    /*Initialize the matrix A that is formulated as a 1-D array*/
    double array[N+2] = {0.0};
    Initiate_Matrix(array);

    /*Problem-1: solve for the maximum and the minimum eigenvalues of the matrix A*/
    double lambda1 = 0.0, lambda501 = 0.0, lambda_s = 0.0;
    double shift = 0.0;
    double eigVector[N] = {1.0};
    //solve for the first eig. value
    Power_algorithm(array, eigVector, &lambda1, shift);
    //printf("The first eig. value obtained is %.11e\n",lambda1);
    //let the shift to be the first lambda_1
    shift = lambda1;
    //solve for the second eig. value
    Power_algorithm(array, eigVector, &lambda501, shift);
    lambda501 = lambda501 + shift;
    //printf("The second eig. value obtained is %.11e\n",lambda501);
    //let lambda1 to be the minimum one while lambda501 the maximum
    if(lambda1>lambda501)
    {
        Exchange(&lambda1,&lambda501);
    }
//    printf("The minimum eigenvalue is: lambda1 = %.11e\n",lambda1);
//    printf("The maximum eigenvalue is: lambda501 = %.11e\n",lambda501);
    //solve for lambda_s by using the inverse Power Alg.
    shift = 0.0;
    inPower_algorithm(array, eigVector, &lambda_s, shift);
//    printf("The min-abs eigenvalue is: lambda_s = %.11e\n",lambda_s);

    /*Problem-2: solve for lambda_i(k), k= 1,2,...,39*/
    int k_num = 39, k;
    double lambda_arr[39] = {1.0};
    for(k=1;k<=k_num;k++)
    {
        shift = lambda1 + k*((lambda501-lambda1)/40.0);
        inPower_algorithm(array, eigVector, &lambda_arr[k-1], shift);
        lambda_arr[k-1] = lambda_arr[k-1] + shift;
//        printf("The eigenvalue lambda_i(%d) is %.11e closed to %.11e\n",k,lambda_arr[k-1],shift);
    }

    /*Problem-3: solve the 2-norm-based condition of the matrix A and its determinant det(A)*/
    double cond = fabs(lambda1/lambda_s);
//    printf("The 2-norm-based condition of this matrix A is cond(A) = %.11e\n", cond);

    double det_value = 0.0;
    shift = 0.0;
    double matrix[R+S+1][N] = {0.0};
    Initiate_CompressedMatrix_shift(matrix, array, shift);
    //step-3: do Doolittle LU decomposition
    Doolittle_LU_composition(matrix, &det_value);//A=LU
//    printf("The determinant of this matrix A is det(A) = %.11e\n", det_value);

    /*record the results*/
    k = 0;
    *(results+k) = lambda1;
    *(results+k+1) = lambda501;
    *(results+k+2) = lambda_s;
    for(k=0;k<39;k++){
        *(results+k+3) = lambda_arr[k];
    }
    *(results+k+3) = cond;
    *(results+k+4) = det_value;

    return 0;
}


/*This function is used to write data into the txt file.*/
int write_data_to_txt(FILE *filepointer, double *array, int n)
{
    int i=0;
    for(i=0;i<n;i++){
        fprintf(filepointer,"%lf\n",*(array+i));
        //printf("Write %lf to the txt file!\n",*(array+i));
    }
    fclose(filepointer);/*close the file pointer*/
    printf("Successfully write the whole data to the txt file!\n");
    return 0;
}


/*test all the sub-functions*/
int test(void)
{
    //initialize the 1-D array which is used to contain all the non-zero elements in A
    double array[N+2] = {0.0};
    Initiate_Matrix(array);

    //test the IndexMapping()
    int i,j;
    for(i=0;i<5;i++){
        for(j=0;j<5;j++){
            printf("A[%d,%d]=%.5f ",i+1,j+1,IndexMapping(array,i,j));
        }
        printf("\n");
    }

    //test the Norm()
    double vector1[N] = {3.,4.};
    printf("The 2-norm of [%.3f,%.3f] is %.3f\n",vector1[0],vector1[1],Norm(vector1));

    //test Multiply_matrix_and_vector()
    double vector2[N] = {1.0};
    Multiply_matrix_and_vector(array,vector1,vector2);
    for(i=0;i<5;i++){
        printf("output[%d]=%.3f\n",i+1,vector2[i]);
    }

    //test Initiate_vector() which is used to randomly initialize each element within [a,b]
    double a=2., b=3.;
    Initiate_vector(vector2, a, b);
    for(i=0;i<5;i++){
        printf("random[%d]=%.3f\n",i+1,vector2[i]);
    }

    //test Initiate_compressedMatrix()
    double matrix[R+S+1][N] = {0.0};
    Initiate_CompressedMatrix(matrix, array);
    for(i=0;i<(R+S+1);i++){
        for(j=0;j<5;j++){
            printf("C[%d,%d]=%.3f ",i+1,j+1,matrix[i][j]);
        }
        printf("\n");
    }

    //test min_two_values() and max_three_values()
    int x=36,y=99,z=88;
    printf("the minimum value among [%d,%d] is %d\n",x,y,min_two_values(x,y));
    printf("the minimum value among [%d,%d,%d] is %d\n",x,y,z,max_three_values(x,y,z));

    //test Doolittle_LU_composition(double **matrix, double *det_value)
    double det_value = 0.0;
    Doolittle_LU_composition(matrix, &det_value);
    printf("The determinant of the matrix A is det(A) = %.11e\n",det_value);

    return 0;
}

/*This is developed based on the Power Algorithm, used to evaluate the
specific abs-max eigenvalue of the given matrix A.
Given that the shift is a double-type variable shift.*/
int Power_algorithm(double *array, double *eigVector, double *eigValue, double shift)
{
    //step-1: initialize a vector u_0 whose elements are uniformly and randomly generated within [minU,maxU]
    double u[N] = {1.0};
    Initiate_vector(u, minU, maxU);
    //given the stop condition flag
    double flag = 1.0, abs_value;
    double temp_eigValue_k_1 = 0.0;
    int i;
    while(flag>epsilon)
    {
        temp_eigValue_k_1 = *(eigValue);
        //calculate ||u(k-1)||
        abs_value = Norm(u);
        if(abs_value<=epsilon){
        printf("The abs-value of u(k) can not be approximatively zero!\n");
        exit(-1);
        }
        //let y(k-1) = u(k-1)/||u(k-1)||
        for(i=0;i<N;i++)
        {
            *(eigVector+i) = (*(u+i))/abs_value;
        }
        //u(k)=Ay(k-1)
        Multiply_matrix_and_vector2(array, eigVector, u, shift);
        *eigValue = innerProduc(eigVector, u);
        flag = sqrt((*eigValue-temp_eigValue_k_1)*(*eigValue-temp_eigValue_k_1))/sqrt((*eigValue)*(*eigValue));
    }
    return 0;
}


/*This is the inverse power algorithm given a shift.*/
int inPower_algorithm(double *array, double *eigVector, double *eigValue, double shift)
{
    //step-1: initialize a vector u_0 whose elements are uniformly and randomly generated within [minU,maxU]
    double u[N] = {1.0};
    double b[N] = {1.0};
    Initiate_vector(u, minU, maxU);
    //step-2: initialize a compressed matrix with a given shift
    double matrix[R+S+1][N] = {0.0};
    Initiate_CompressedMatrix_shift(matrix, array, shift);
    //step-3: do Doolittle LU decomposition
    double det_value = 0.0;
    Doolittle_LU_composition(matrix, &det_value);//A=LU
    //given the stop condition flag
    double flag = 1.0, abs_value;
    double temp_eigValue_k_1 = 0.0;
    int i;
    while(flag>epsilon)
    {
        temp_eigValue_k_1 = *(eigValue);
        //||u(k-1)||
        abs_value = Norm(u);
        if(abs_value<=epsilon){
        printf("The abs-value of u(k) can not be approximatively zero!\n");
        exit(-1);
        }
        //let y(k-1) = u(k-1)/||u(k-1)||
        for(i=0;i<N;i++)
        {
            *(eigVector+i) = (*(u+i))/abs_value;
            *(b+i) = *(eigVector+i);
        }
        //solve Au(k)=y(k-1) to get u(k)
        Solve_linear_equation(matrix,b,u);
        *(eigValue) = 1.0/innerProduc(eigVector, u);
        flag = sqrt((*eigValue-temp_eigValue_k_1)*(*eigValue-temp_eigValue_k_1))/sqrt((*eigValue)*(*eigValue));
    }

    return 0;
}


/*This function calculates the product of the matrix A and the given vector.*/
int Multiply_matrix_and_vector(double *array, double *vector, double *output)
{
    int i,j;
    double temp;
    for(i=0;i<N;i++){
        temp = 0.0;
        for(j=0;j<N;j++){
            temp += (IndexMapping(array,i,j))*(*(vector+j));
        }
        *(output+i) = temp;
    }
    return 0;
}


/*This function calculates the product of the matrix A and the given vector given a shift.*/
int Multiply_matrix_and_vector2(double *array, double *vector, double *output, double shift)
{
    int i,j;
    double temp;
    for(i=0;i<N;i++){
        temp = 0.0;
        for(j=0;j<N;j++){
            if(i==j)
            {
                temp += (IndexMapping(array,i,j)-shift)*(*(vector+j));//a_ii - p
            }
            else
            {
                temp += (IndexMapping(array,i,j))*(*(vector+j));//a_ij
            }
        }
        *(output+i) = temp;
    }
    return 0;
}


/*This function calculates the inner-product of two given vectors.*/
double innerProduc(double *arr1, double *arr2)
{
    double temp = 0.0;
    int i;
    for(i=0;i<N;i++){
        temp += (*(arr1+i))*(*(arr2+i));
    }
    return temp;
}


/*This function calculates the 2-norm of a vector*/
double Norm(double *array)
{
    double temp=0.0;
    int i;
    for(i=0;i<N;i++){
        temp += (*(array+i))*(*(array+i));
    }
    return sqrt(temp);
}


/*This function is used to specify the real symmetric matrix A, in which
1-D array is used to store all the different non-zero elements in A in order to
reduce the memory occupancy. Note that the number of different non-zero elements
is (N+2).
Given that
arr is a double-type pointer associate with an array whose length is set to (N+2)*/
int Initiate_Matrix(double *arr)
{
    int k;
    for(k=0;k<N;k++)
    {
        *(arr+k) = (1.64-0.024*(k+1))*sin(0.2*(k+1))-0.64*exp((0.1)/(1.0*(k+1)));//a_i
    }
    *(arr+k) = 0.16;
    *(arr+k+1) = -0.064;
    return 0;
}


/*This function initiates a vector by randomness.*/
int Initiate_vector(double *array, double a, double b)
{
    srand( (unsigned)time( NULL ) );
    int i;
    for(i=0;i<N;i++){
        *(array+i) = 1.0*a+(b-a)*(rand()/32767.0);
    }
    return 0;
}

/*This function is used to initiate a compressed matrix applied to store
the non-zero elements in the original matrix A. Note that the non-zero
elements are mainly located within the bandwidth of the matrix A.
Given that the compressed array is a 1-D array containing all the non-zero
elements of A.*/
int Initiate_CompressedMatrix(double **matrix, double *array)
{
    int i, j;
    int ci, cj;
    int cn=N;//the column number of the compressed matrix, cn = N
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(((i+S)>=j)&&(i<=(j+R)))//the A(i,j) is within the non-zero bandwidth
            {
                ci = (i-j+S+1)-1;
                cj = j;
                *((double*)matrix+ci*cn+cj)=IndexMapping(array, i, j);//C(i-j+S+1,j)=A(i,j)
            }
        }
    }
    return 0;
}


/*This function is used to initiate a compressed matrix applied to store
the non-zero elements in the original matrix A. Note that the non-zero
elements are mainly located within the bandwidth of the matrix A.
Given that the compressed array is a 1-D array containing all the non-zero
elements of A. Given a shift.*/
int Initiate_CompressedMatrix_shift(double **matrix, double *array, double shift)
{
    int i, j;
    int ci, cj;
    int cn=N;//the column number of the compressed matrix, cn = N
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(((i+S)>=j)&&(i<=(j+R)))//the A(i,j) is within the non-zero bandwidth
            {
                ci = (i-j+S+1)-1;
                cj = j;
                if (i==j)
                {
                    *((double*)matrix+ci*cn+cj)=IndexMapping(array, i, j)-shift;//C(i-j+S+1,j)=A(i,j)
                }
                else
                {
                    *((double*)matrix+ci*cn+cj)=IndexMapping(array, i, j);//C(i-j+S+1,j)=A(i,j)
                }

            }
        }
    }
    return 0;
}


/*This function is used to map the 2-D index (i,j) of the matrix A to the corresponding index
in the array arr. With this function, we can get the element A(i,j) when given arr, i, and j.
Note that in C, the initial index is 0, so that the range of i or j should be [0,N-1].*/
double IndexMapping(double *arr, int i, int j)
{
    double element = 0.0;
    //int bandwidth = (2*2 + 1);//the bandwidth of the sparse matrix A
    if(i == j)
    {
        element = *(arr+i);//a_i
    }
    else if(((i+1) == j)||(i == (j+1)))
    {
        element = *(arr+N+2-2);//b
    }
    else if(((i+2) == j)||(i == (j+2)))
    {
        element = *(arr+N+2-1);//c
    }
    return element;
}


/*This function exchanges the two values of the given parameters using their addresses.*/
int Exchange(double *Value1, double *Value2)
{
    double temp;
    temp = *Value1;
    *Value1 = *Value2;
    *Value2 = temp;
    return 0;
}


/*This function determine which one is the maximum/minimum value given
any two double variables.*/
int Determine_Max_Min(double *maxV, double *minV)
{
    if((*maxV)<(*minV))
    {
        Exchange(maxV, minV);
    }
    return 0;
}


/*This function solve the linear equation Au(k)=y(k-1) based on the Doolittle LU decomposition*/
int Solve_linear_equation(double **matrix, double *b, double *solution)
{
    double det_value = 0.0;

    int t_lower, t_upper, i;

    for(i=2;i<=N;i++)
    {
        t_lower = (i-R)>1?(i-R):1;
        t_upper = i-1;
        *(b+i-1) = *(b+i-1) - Sum_up3(matrix,b,i,t_lower,t_upper);
    }
    *(solution + N-1) = *(b+N-1)/(*((double*)matrix+N*S+N-1));
    for(i=N-1;i>=1;i--)
    {
        t_lower = i+1;
        t_upper = (i+S)>N?N:(i+S);
        *(solution+i-1) = ((*(b+i-1))-Sum_up4(matrix,solution,i,t_lower,t_upper))/(*((double*)matrix+N*S+i-1));
    }
    return 0;
}


/*This function sums up all the c(i-t+s+1,t)x(t) with respect to t.*/
double Sum_up4(double **matrix, double *solution, int i, int t_lower, int t_upper)
{
    double temp=0.0;
    int ci, cj, t;
    for(t=t_lower;t<=t_upper;t++)
    {
        ci = (i-t+S+1)-1;
        cj = t-1;
        temp += (*((double*)matrix+N*ci+cj))*(*(solution+t-1));
    }
    return temp;
}


/*This function sums up all the c(i-t+s+1,t)b(t) with respect to t.*/
double Sum_up3(double **matrix,double *b,int i,int t_lower,int t_upper)
{
    double temp=0.0;
    int ci, cj, t;
    for(t=t_lower;t<=t_upper;t++)
    {
        ci = (i-t+S+1)-1;
        cj = t-1;
        temp += (*((double*)matrix+N*ci+cj))*(*(b+t-1));
    }
    return temp;
}


/*This function calculates the determinant of the matrix based on
the Doolittle LU composition algorithm:
the following function is developed based on
the well-know Doolittle LU composition algorithm
to solve a linear equations system. It should be
noted that this is a basic Doolittle algorithm, and it
requires that the diagonal elements cannot be zero.*/
int Doolittle_LU_composition(double **matrix, double *det_value)
{
    int k, i, j, t, t_lower, t_upper, i_upper, j_upper;
    int ci, cj;
    double temp = 1.0;
    for(k=1;k<=N;k++)
    {
        //calculate the U matrix
        j_upper=min_two_values(k+S,N);
        for(j=k;j<=j_upper;j++)
        {
            t_lower = max_three_values(1,k-R,j-S);
            t_upper = k-1;
            ci = (k-j+S+1)-1;
            cj = j-1;
            *((double*)matrix+N*ci+cj) = *((double*)matrix+N*ci+cj) - Sum_up(matrix, t_lower, t_upper, k, j);

        }

        //calculate the L matrix
        if(k<N){
            i_upper = min_two_values(k+R,N);
            for(i=(k+1);i<=i_upper;i++){
                t_lower = max_three_values(1,i-R,k-S);
                t_upper = k-1;
                ci = (i-k+S+1)-1;
                cj = k-1;
                *((double*)matrix+N*ci+cj) = (*((double*)matrix+N*ci+cj)-Sum_up2(matrix, t_lower, t_upper, k, i))/(*((double*)matrix+N*S+k-1));
            }
        }
    }
    for(k=1;k<=N;k++)
    {
        ci = S;
        cj = k - 1;
        temp *= *((double*)matrix+N*ci+cj);
    }
    *det_value = temp;
    return 0;

}


/*This function sums up the product of c(k-t+s+1,t)*c(t-j+s+1,j)
given the lower bound of t, max{1,k-r,j-s}, and the upper bound k-1*/
double Sum_up(double **matrix, int t1, int t2, int k, int j)
{
    double temp = 0.0;
    int t, ci1, cj1, ci2, cj2;
    for(t=t1;t<=t2;t++)
    {
        ci1 = (k-t+S+1)-1;
        cj1 = t-1;
        ci2 = (t-j+S+1)-1;
        cj2 = j-1;
        temp += (*((double*)matrix+N*ci1+cj1))*(*((double*)matrix+N*ci2+cj2));
    }
    return temp;
}


/*This function sums up the inner elements.*/
double Sum_up2(double ** matrix, int t1, int t2, int k, int i)
{
    double temp = 0.0;
    int t, ci1, cj1, ci2, cj2;
    for(t=t1;t<=t2;t++){
        ci1 = (i-t+S+1)-1;
        cj1 = t-1;
        ci2 = (t-k+S+1)-1;
        cj2 = k-1;
        temp += (*((double*)matrix+N*ci1+cj1))*(*((double*)matrix+N*ci2+cj2));
    }
    return temp;
}


/*This function selects the minimum one from two given variables.*/
int min_two_values(int a, int b)
{
    return a<b?a:b;
}


/*This function selects the maximum value from three given variables.*/
int max_three_values(int a, int b, int c)
{
    int temp = a>=b?a:b;
    return temp>c?temp:c;
}
