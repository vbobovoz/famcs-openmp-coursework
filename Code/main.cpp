#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "omp.h"
#include <cassert>
#include <chrono>

void solveTridiagonalSystem(int i,
                            const std::vector<std::vector<double>> &cC,
                            const std::vector<std::vector<double>> &aC,
                            const std::vector<std::vector<double>> &bC,
                            std::vector<double> &phi,
                            std::vector<double> &Y_i,
                            int N_2)
{
    int m = N_2 - 1;
    std::vector<double> alpha(m+1,0.0);
    std::vector<double> beta(m+1,0.0);

    // Инициализация первых коэффициентов прогонки
    double gamma = cC[i][0];
    alpha[1] = -bC[i][0]/gamma;
    beta[1]  = phi[0]/gamma;

    // Прямая прогонка для j=1..m-2
    for (int j=1; j<=m-2; j++) {
        gamma = cC[i][j] + aC[i][j]*alpha[j];
        alpha[j+1] = -bC[i][j]/gamma;
        beta[j+1]  = (phi[j]-aC[i][j]*beta[j])/gamma;
    }

    // Последний элемент
    gamma = cC[i][m-1] + aC[i][m-1]*alpha[m-1];
    beta[m] = (phi[m-1]-aC[i][m-1]*beta[m-1])/gamma;

    // Обратная прогонка
    Y_i[m-1] = beta[m];
    for (int j=m-2; j>=0; j--) {
        Y_i[j] = alpha[j+1]*Y_i[j+1] + beta[j+1];
    }
}

void initializeSystem(int N_1, int N_2,
                      std::vector<std::vector<double>> &cC,
                      std::vector<std::vector<double>> &aC,
                      std::vector<std::vector<double>> &bC,
                      std::vector<std::vector<double>> &cA,
                      std::vector<std::vector<double>> &aA,
                      std::vector<std::vector<double>> &bA,
                      std::vector<std::vector<double>> &cB,
                      std::vector<std::vector<double>> &aB,
                      std::vector<std::vector<double>> &bB,
                      std::vector<std::vector<double>> &fslae,
                      std::vector<double> &F)
{
    // Размер для внутренних массивов: (N_1-1)x(N_2-1)
    cC.assign(N_1-1,std::vector<double>(N_2-1,0.0));
    aC.assign(N_1-1,std::vector<double>(N_2-1,0.0));
    bC.assign(N_1-1,std::vector<double>(N_2-1,0.0));

    cA.assign(N_1-1,std::vector<double>(N_2-1,0.0));
    aA.assign(N_1-1,std::vector<double>(N_2-1,0.0));
    bA.assign(N_1-1,std::vector<double>(N_2-1,0.0));

    cB.assign(N_1-1,std::vector<double>(N_2-1,0.0));
    aB.assign(N_1-1,std::vector<double>(N_2-1,0.0));
    bB.assign(N_1-1,std::vector<double>(N_2-1,0.0));

    fslae.assign(N_1-1,std::vector<double>(N_2-1,0.0));

    double denom = 2.0*N_1 + N_2;
    for (int i=1; i<=N_1-1; i++) {
        for (int j=1; j<=N_2-1; j++) {
            // C_i
            cC[i-1][j-1] = 8.0;
            if (j>=2) aC[i-1][j-1]=(2.0*i+j)/denom;
            if (j<=N_2-2) bC[i-1][j-1]=(2.0*i+j)/denom;

            // A_i
            cA[i-1][j-1]=(2.0*i+j)/denom;
            if (j>=2) aA[i-1][j-1]=(2.0*i+j)/denom;
            if (j<=N_2-2) bA[i-1][j-1]=(2.0*i+j)/denom;

            // B_i
            cB[i-1][j-1]=(2.0*i+j)/denom;
            if (j>=2) aB[i-1][j-1]=(2.0*i+j)/denom;
            if (j<=N_2-2) bB[i-1][j-1]=(2.0*i+j)/denom;
        }
    }

    // Формируем F для решения (1,...,1)
    int N = (N_1+1)*(N_2-1);
    F.assign(N,0.0);

    // F_0=(1,...,1)
    for (int j=1; j<=N_2-1; j++) {
        F[(0)*(N_2-1)+(j-1)] = 1.0;
    }

    // F_(N_1)=(1,...,1)
    for (int j=1; j<=N_2-1; j++) {
        F[(N_1)*(N_2-1)+(j-1)] = 1.0;
    }

    // Внутренние F_i
    for (int i=1; i<=N_1-1; i++) {
        for (int j=1; j<=N_2-1; j++) {
            double val=0.0;
            // A_i
            val+= cA[i-1][j-1];
            if (j>=2) val+=aA[i-1][j-1];
            if (j<=N_2-2) val+=bA[i-1][j-1];

            // C_i
            val+= cC[i-1][j-1];
            if (j>=2) val+=aC[i-1][j-1];
            if (j<=N_2-2) val+=bC[i-1][j-1];

            // B_i
            val+= cB[i-1][j-1];
            if (j>=2) val+=aB[i-1][j-1];
            if (j<=N_2-2) val+=bB[i-1][j-1];

            // F_i
            F[i*(N_2-1)+(j-1)] = val;
        }
    }

    // fslae(i,j) = F_i(j) для i=1..N_1-1
    for (int i=1; i<=N_1-1; i++) {
        for (int j=1; j<=N_2-1; j++) {
            fslae[i-1][j-1] = F[i*(N_2-1)+(j-1)];
        }
    }
}

void gaussSeidel(int N_1, int N_2, int K,
                 const std::vector<std::vector<double>> &cC,
                 const std::vector<std::vector<double>> &aC,
                 const std::vector<std::vector<double>> &bC,
                 const std::vector<std::vector<double>> &cA,
                 const std::vector<std::vector<double>> &aA,
                 const std::vector<std::vector<double>> &bA,
                 const std::vector<std::vector<double>> &cB,
                 const std::vector<std::vector<double>> &aB,
                 const std::vector<std::vector<double>> &bB,
                 const std::vector<std::vector<double>> &fslae,
                 const std::vector<double> &F,
                 std::vector<std::vector<double>> &y)
{
    // Начальные условия
    // Y_0 = F_0, Y_(N_1)=F_(N_1)
    for (int j=0; j<N_2-1; j++) {
        y[0][j] = F[j];
        y[N_1][j] = F[N_1*(N_2-1)+j];
    }

    // Внутренние y(i,j)=0
    for (int i=1; i<=N_1-1; i++) {
        for (int j=0; j<N_2-1; j++)
            y[i][j]=0.0;
    }

    std::vector<double> phi(N_2-1,0.0);
    std::vector<double> Y_i(N_2-1,0.0);

    for (int k_iter=0; k_iter<K; k_iter++) {
        for (int i=1; i<=N_1-1; i++) {
            for (int j=1; j<=N_2-1; j++) {
                double val = fslae[i-1][j-1];

                int jidx=j-1;
                double y_im1_j   = y[i-1][jidx];
                double y_im1_jm1 = (j>1)     ? y[i-1][jidx-1]:0.0;
                double y_im1_jp1 = (j<N_2-1) ? y[i-1][jidx+1]:0.0;

                double y_ip1_j   = y[i+1][jidx];
                double y_ip1_jm1 = (j>1)     ? y[i+1][jidx-1]:0.0;
                double y_ip1_jp1 = (j<N_2-1) ? y[i+1][jidx+1]:0.0;

                // Вычитаем A_i Y_(i-1)
                val -= cA[i-1][j-1]*y_im1_j;
                if (j>=2) val -= aA[i-1][j-1]*y_im1_jm1;
                if (j<=N_2-2) val -= bA[i-1][j-1]*y_im1_jp1;

                // Вычитаем B_i Y_(i+1)
                val -= cB[i-1][j-1]*y_ip1_j;
                if (j>=2) val -= aB[i-1][j-1]*y_ip1_jm1;
                if (j<=N_2-2) val -= bB[i-1][j-1]*y_ip1_jp1;

                phi[j-1]=val;
            }

            solveTridiagonalSystem(i-1,cC,aC,bC,phi,Y_i,N_2);

            for (int j=0; j<N_2-1; j++) {
                y[i][j]=Y_i[j];
            }
        }
    }
}

// Максимум-норма ошибки: max|x_k - 1|
double maxNormError(const std::vector<double> &x) {
    double max_error = 0.0;
    for (double val : x) {
        double err = std::fabs(val - 1.0);
        if (err > max_error)
            max_error = err;
    }
    return max_error;
}

void printMatrixBlock(const std::vector<std::vector<double>> &mat) {
    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            std::cout << std::setw(10) << mat[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void printVector(const std::vector<double> &vec) {
    for (double val : vec) {
        std::cout << val << " ";
    }
    std::cout << "\n";
}

void putBlock(std::vector<std::vector<double>> &FullMat,
              int N_1, int N_2,
              int blockRow, int blockCol,
              const std::vector<std::vector<double>> *cM,
              const std::vector<std::vector<double>> *aM,
              const std::vector<std::vector<double>> *bM,
              bool isIdentity = false,
              bool isZero = false)
{
    int rowStart = blockRow*(N_2-1);
    int colStart = blockCol*(N_2-1);

    int m = N_2-1; 
    if (isIdentity) {
        for (int j=0; j<m; j++) {
            FullMat[rowStart+j][colStart+j] = 1.0;
        }
        return;
    }
    if (isZero) {
        // Нулевой блок - ничего не делаем
        return;
    }

    // обычный трехдиагональный блок
    // i_index = blockRow-1 для внутренних блоков (A_i, B_i, C_i)
    // но если blockRow=0 или blockRow=N_1, значит это граничные блоки.
    int i_index = blockRow-1; 
    for (int j=1; j<=m; j++) {
        int rr = rowStart + (j-1);
        double c_val = (*cM)[i_index][j-1];
        FullMat[rr][colStart+(j-1)] = c_val;

        if (j>=2) {
            double a_val = (*aM)[i_index][j-1];
            FullMat[rr][colStart+(j-1)-1] = a_val;
        }
        if (j<=m-1) {
            double b_val = (*bM)[i_index][j-1];
            FullMat[rr][colStart+(j-1)+1] = b_val;
        }
    }
}

void printFullMatrix(int N_1, int N_2,
                     const std::vector<std::vector<double>> &cC,
                     const std::vector<std::vector<double>> &aC,
                     const std::vector<std::vector<double>> &bC,
                     const std::vector<std::vector<double>> &cA,
                     const std::vector<std::vector<double>> &aA,
                     const std::vector<std::vector<double>> &bA,
                     const std::vector<std::vector<double>> &cB,
                     const std::vector<std::vector<double>> &aB,
                     const std::vector<std::vector<double>> &bB)
{
    int N = (N_1+1)*(N_2-1);
    std::vector<std::vector<double>> FullMat(N, std::vector<double>(N,0.0));

    // Первая строка блоков: C_0=I, B_0=0
    putBlock(FullMat, N_1,N_2, 0,0,nullptr,nullptr,nullptr,true,false); // C_0=I
    if (N_1>=1) putBlock(FullMat, N_1,N_2, 0,1,nullptr,nullptr,nullptr,false,true); // B_0=0

    // i=1..N_1-1: [A_i, C_i, B_i]
    for (int i=1; i<=N_1-1; i++) {
        putBlock(FullMat,N_1,N_2, i,i-1,&cA,&aA,&bA);
        putBlock(FullMat,N_1,N_2, i,i,  &cC,&aC,&bC);
        putBlock(FullMat,N_1,N_2, i,i+1,&cB,&aB,&bB);
    }

    // i=N_1: A_(N_1)=0, C_(N_1)=I
    if (N_1>=1) {
        putBlock(FullMat,N_1,N_2, N_1, N_1-1,nullptr,nullptr,nullptr,false,true); // A_(N_1)=0
        putBlock(FullMat,N_1,N_2, N_1, N_1,   nullptr,nullptr,nullptr,true,false); // C_(N_1)=I
    }

    // Выводим полную матрицу
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            std::cout << std::setw(10) << FullMat[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void gaussSeidelParallel(int N_1, int N_2, int K,
                         const std::vector<std::vector<double>> &cC,
                         const std::vector<std::vector<double>> &aC,
                         const std::vector<std::vector<double>> &bC,
                         const std::vector<std::vector<double>> &cA,
                         const std::vector<std::vector<double>> &aA,
                         const std::vector<std::vector<double>> &bA,
                         const std::vector<std::vector<double>> &cB,
                         const std::vector<std::vector<double>> &aB,
                         const std::vector<std::vector<double>> &bB,
                         const std::vector<std::vector<double>> &fslae,
                         const std::vector<double> &F,
                         std::vector<std::vector<double>> &y)
{
    // Инициализация
    for (int j=0;j<N_2-1;j++){
        y[0][j]=F[j];
        y[N_1][j]=F[N_1*(N_2-1)+j];
    }
    for (int i=1;i<=N_1-1;i++){
        for (int j=0;j<N_2-1;j++)
            y[i][j]=0.0;
    }

    int m=N_2-1;
    std::vector<std::vector<double>> Y_temp(N_1+1,std::vector<double>(m,0.0));

    for (int k_iter=0;k_iter<K;k_iter++){
        std::vector<std::vector<double>> Y_old=y;

        // i,j не указываем в private, т.к. они будут объявлены как итераторы цикла for,
        // что делает их автоматически private для данного цикла.
        #pragma omp parallel default(none) \
            shared(N_1,N_2,m,Y_temp,y,Y_old,fslae,F,cC,aC,bC,cA,aA,bA,cB,aB,bB)
        {
            // Локальные для каждого потока переменные
            std::vector<double> phi(m,0.0);
            std::vector<double> Y_i(m,0.0);
            double val;
            int jidx;

            #pragma omp for
            for (int i=1;i<=N_1-1;i++){
                for (int j=1;j<=m;j++){
                    val=fslae[i-1][j-1];

                    jidx=j-1;
                    double y_im1_j=Y_old[i-1][jidx];
                    double y_im1_jm1=(j>1)?Y_old[i-1][jidx-1]:0.0;
                    double y_im1_jp1=(j<m)?Y_old[i-1][jidx+1]:0.0;

                    double y_ip1_j=Y_old[i+1][jidx];
                    double y_ip1_jm1=(j>1)?Y_old[i+1][jidx-1]:0.0;
                    double y_ip1_jp1=(j<m)?Y_old[i+1][jidx+1]:0.0;

                    val -= cA[i-1][j-1]*y_im1_j;
                    if (j>=2) val -= aA[i-1][j-1]*y_im1_jm1;
                    if (j<=m-1) val -= bA[i-1][j-1]*y_im1_jp1;

                    val -= cB[i-1][j-1]*y_ip1_j;
                    if (j>=2) val -= aB[i-1][j-1]*y_ip1_jm1;
                    if (j<=m-1) val -= bB[i-1][j-1]*y_ip1_jp1;

                    phi[j-1]=val;
                }

                solveTridiagonalSystem(i-1,cC,aC,bC,phi,Y_i,N_2);

                for (int jj=0;jj<m;jj++){
                    Y_temp[i][jj]=Y_i[jj];
                }
            }
        }

        for (int i=1;i<=N_1-1;i++){
            for (int j=0;j<m;j++){
                y[i][j]=Y_temp[i][j];
            }
        }
    }
}

int main() {
    int N_1=999; // N_1 + 1 = количество блоков в строке и столбце
    int N_2=1001; // N_2 - 1 = размерность блока (количество строк и столбцов в блоке)
    int K=2; // Число итераций

    std::vector<std::vector<double>> cC1,aC1,bC1,cA1,aA1,bA1,cB1,aB1,bB1,fslae1;
    std::vector<double> F1;

    // Инициализация системы
    initializeSystem(N_1,N_2,cC1,aC1,bC1,cA1,aA1,bA1,cB1,aB1,bB1,fslae1,F1);

    // Выводим вектор правой части
    // std::cout << "Right-hand side vector F:\n";
    // printVector(F1);
    // std::cout << "\n";

    // Выводим один из блоков
    // std::cout << "\nBlock aC:\n";
    // printMatrixBlock(aC1);
    // std::cout << "\n";

    // Вывод полной матрицы
    // std::cout << "\nFull Matrix:\n";
    // printFullMatrix(N_1,N_2,cC1,aC1,bC1,cA1,aA1,bA1,cB1,aB1,bB1);
    // std::cout << "\n";

    // Массив решения y(i,j)
    std::vector<std::vector<double>> y1(N_1+1, std::vector<double>(N_2-1,0.0));

    // Запуск метода Гаусса–Зейделя
    auto start_seq=std::chrono::high_resolution_clock::now();
    gaussSeidel(N_1,N_2,K,cC1,aC1,bC1,cA1,aA1,bA1,cB1,aB1,bB1,fslae1,F1,y1);
    auto end_seq=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> seq_dur=end_seq-start_seq;
    double seq_time1=seq_dur.count();
    
    std::cout << "Sequential time: " << seq_time1 << " s\n";
    
    // Двумерный массив решения y в одномерный вектор x
    std::vector<double> x1((N_1+1)*(N_2-1),0.0);
    for (int i=0; i<=N_1; i++) {
        for (int j=0; j<N_2-1; j++) {
            int idx = i*(N_2-1) + j;
            x1[idx] = y1[i][j];
        }
    }

    // std::cout << "Solution vector-1:\n";
    // printVector(x1);

    // Вычисляем максимум-норму ошибки
    double error1 = maxNormError(x1);
    std::cout << "Maximum norm of error: " << error1 << "\n\n";

    // -------------- Параллельный метод Гаусса-Зейделя -------------- //
    std::vector<std::vector<double>> cC2,aC2,bC2,cA2,aA2,bA2,cB2,aB2,bB2,fslae2;
    std::vector<double> F2;

    // Инициализация системы
    initializeSystem(N_1,N_2,cC2,aC2,bC2,cA2,aA2,bA2,cB2,aB2,bB2,fslae2,F2);

    // Выводим вектор правой части
    // std::cout << "Right-hand side vector F:\n";
    // printVector(F2);
    // std::cout << "\n";

    std::vector<std::vector<double>> y2(N_1+1,std::vector<double>(N_2-1,0.0));

    start_seq=std::chrono::high_resolution_clock::now();
    gaussSeidelParallel(N_1,N_2,K,cC2,aC2,bC2,cA2,aA2,bA2,cB2,aB2,bB2,fslae2,F2,y2);
    end_seq=std::chrono::high_resolution_clock::now();
    seq_dur=end_seq-start_seq;
    double seq_time2=seq_dur.count();
    std::cout << "Parallel time: " << seq_time2 << " s\n";
    
    // Двумерный массив решения y в одномерный вектор x
    std::vector<double> x2((N_1+1)*(N_2-1),0.0);
    for (int i=0; i<=N_1; i++) {
        for (int j=0; j<N_2-1; j++) {
            int idx = i*(N_2-1) + j;
            x2[idx] = y2[i][j];
        }
    }

    // std::cout << "Solution vector-2:\n";
    // printVector(x2);

    // Вычисляем максимум-норму ошибки
    double error2 = maxNormError(x2);
    std::cout << "Maximum norm of error: " << error2 << "\n\n";

    // Ускорение
    std::cout << "Acceleration: " << seq_time1/seq_time2 << "\n";

    return 0;
}