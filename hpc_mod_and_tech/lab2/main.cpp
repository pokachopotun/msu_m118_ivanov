#define _USE_MATH_DEFINES
#include <cstdio>

#include <cmath>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>

using namespace std;

double u1(double x, double y, double t) {
    return (t * t  + 1) * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

double u1t(double x, double y, double t) {
    return 2 * t * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

double u1xx(double x, double y, double t) {
    return (t * t  + 1) * sin(2 * M_PI * x) * sin(2 * M_PI * y) * (-4 * M_PI * M_PI);
}

double u1yy(double x, double y, double t) {
    return (t * t  + 1) * sin(2 * M_PI * x) * sin(2 * M_PI * y) * (-4 * M_PI * M_PI);
}

double u1xy(double x, double y, double t) {
    return (t * t  + 1) * cos(2 * M_PI * x) * cos(2 * M_PI * y) * (4 * M_PI * M_PI);
}
        
double u2(double x, double y, double t) {
    return (t * t + 1) * (cos(2 * M_PI * x) - 1) * (cos(2 * M_PI * y) - 1);
}

double u2t(double x, double y, double t) {
    return 2 * t * (cos(2 * M_PI * x) - 1) * (cos(2 * M_PI * y) - 1);
}

double u2xx(double x, double y, double t) {
    return (t * t + 1) * cos(2 * M_PI * x) * (cos(2 * M_PI * y) - 1) * (-4 * M_PI * M_PI);
}

double u2yy(double x, double y, double t) {
    return (t * t + 1) * cos(2 * M_PI * y) * (cos(2 * M_PI * x) - 1) * (-4 * M_PI * M_PI);
}

double u2xy(double x, double y, double t) {
    return (t * t + 1) * sin(2 * M_PI * x) * sin(2 * M_PI * y) * (4 * M_PI * M_PI);
}

double f1(double x, double y, double t) {
    return u1t(x,y,t) - 2 * u1xx(x,y,t) - u1yy(x,y,t) - u2xy(x,y,t);
}

double f2(double x, double y, double t) {
    return u2t(x,y,t) - u2xx(x,y,t) - 2 * u2yy(x,y,t) - u1xy(x,y,t);
}

class Calc {
private:
    int tl, t, tr;
    int l, m, r;
    int bl, b, br;
    double hx, hy, ht;
    int i, j;
    double x, y, time;

    double Laplas(double* v) {
        return (v[t] - 2 * v[m] + v[b]) / (hx * hx) + (v[l] - 2 * v[m] + v[r]) / (hy * hy);
    }

    double Grad1(double* v1, double* v2) {
        return (v1[t] - 2 * v1[m] + v1[b]) / (hx * hx) + (v2[br] - v2[bl] - v2[tr] + v2[tl]) / (4 * hx * hy);
    }

    double Grad2(double* v1, double* v2) {
        return (v2[l] - 2 * v2[m] + v2[r]) / (hy * hy) + (v1[br] - v1[bl] - v1[tr] + v1[tl]) / (4 * hx * hy);
    }

public:
    Calc(int localYSize, double hx, double hy, double ht,  int i, int j, double x, double y, double time) 
            : hx(hx), hy(hy), ht(ht), i(i), j(j), x(x), y(y), time(time)
    {
        int storedYSize = localYSize + 2; 

        tl = storedYSize * (i) + j;
        t = storedYSize * (i) + 1 + j;
        tr = storedYSize * (i) + 1 + j + 1;

        l = storedYSize * (i + 1) + j;
        m = storedYSize * (i + 1) + 1 + j;
        r = storedYSize * (i + 1) + 1 + j + 1;

        bl = storedYSize * (i + 2) + j;
        b = storedYSize * (i + 2) + 1 + j;
        br = storedYSize * (i + 2) + 1 + j + 1;
    }

    void CalcNextV1(double* next1, double* v1, double* v2) {
        double laplas1 = Laplas(v1);      
        double g1 = Grad1(v1, v2);
        next1[m] = ht * (g1 + laplas1 + f1(x,y,time)) + v1[m];
    }

    void CalcNextV2(double* next2, double* v1, double* v2) {
        double laplas2 = Laplas(v2);
        double g2 = Grad2(v1, v2);
        next2[m] = ht * (g2 + laplas2 + f2(x,y,time)) + v2[m];
    }
};

class Solution {
private:
    // mpi grid size
    int mpiN, mpiM;
    // solution grid size 
    int gridN, gridM, gridT;
    // constraints
    double xmax, ymax;
    // grid step
    double hx, hy, ht;
    // mpi params
    int mpiNodes, mpiRank;
    int mpiX, mpiY, isXMax, isYMax;
    int localXSize, localYSize, storedXSize, storedYSize;
    int sliceSize;
    int commBufSize;
    int xStart, yStart;
    vector<double*> v1;
    vector<double*> v2;
    double* sendBuf;
    double* recvBuf;

    int Neighbor[8];
    int Tag[8];
    int CommBufShift[8];
    int CommBufSize[8];

    int cyclicBufSize;
    double* cyclicSendBuf;
    double* cyclicRecvBuf;

    int ompNumThreads; 
    int printDebug;

    void Init(int argc, char** argv) {
        if (argc < 9) {
            cerr << "use ./solution mpiN mpiM gridN gridM gridT xmax ymax ompNumThreads printDebug" << endl;
            return;
        }
        // mpi grid size
        mpiN = atoi(argv[1]);
        mpiM = atoi(argv[2]);

        // solution grid size 
        gridN = atoi(argv[3]);
        gridM = atoi(argv[4]);
        gridT = atoi(argv[5]);
        gridT += 1;

        // constraints
        xmax = atof(argv[6]);
        ymax = atof(argv[7]);

        ompNumThreads = atoi(argv[8]);

        omp_set_num_threads(ompNumThreads);

        printDebug = 0;
        if (argc >= 10) {
            printDebug = atoi(argv[9]);
        }

        // grid step
        hx = xmax / (gridN - 1);
        hy = ymax / (gridM - 1);
        ht = min(hx * hx, hy * hy) / 1000; 

        MPI_Comm_size(MPI_COMM_WORLD, &mpiNodes);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        if (mpiNodes != mpiN * mpiM) {
            cerr << "mpiNodes != mpiN * mpiM" << endl;
            return;
        }

        mpiX = mpiRank / mpiM;
        mpiY = mpiRank % mpiM;

        bool isXMax = mpiX == mpiN - 1;  
        bool isYMax = mpiY == mpiM - 1;

        int maxLocalXSize = gridN / mpiN;
        localXSize = isXMax ? maxLocalXSize + gridN % maxLocalXSize : maxLocalXSize;

        int maxLocalYSize = gridM / mpiM;
        localYSize = isYMax ? maxLocalYSize + gridM % maxLocalYSize : maxLocalYSize;

        storedXSize = localXSize + 2;
        storedYSize = localYSize + 2;

        sliceSize = storedXSize * storedYSize;
        commBufSize = 2 * (localXSize + localYSize) + 4;

        xStart = mpiX * maxLocalXSize;
        yStart = mpiY * maxLocalYSize;

        v1.reserve(2);
        v2.reserve(2);
        
        for (int i = 0; i < 2; i++) {
            v1.push_back(new double[sliceSize]);
            v2.push_back(new double[sliceSize]);
        }

        sendBuf = new double[commBufSize];
        recvBuf = new double[commBufSize];

        cyclicBufSize = 2 * localXSize;
        cyclicSendBuf = new double[cyclicBufSize];
        cyclicRecvBuf = new double[cyclicBufSize];

        Neighbor[0] = GetDest(-1, 0);
        Neighbor[1] = GetDest(1, 0);
        Neighbor[2] = GetDest(0, -1);
        Neighbor[3] = GetDest(0, 1);
        Neighbor[4] = GetDest(-1, -1);
        Neighbor[5] = GetDest(-1, 1);
        Neighbor[6] = GetDest(1, -1);
        Neighbor[7] = GetDest(1, 1);

        Tag[0] = 1;
        Tag[1] = 0;
        Tag[2] = 3;
        Tag[3] = 2;
        Tag[4] = 7;
        Tag[5] = 6;
        Tag[6] = 5;
        Tag[7] = 4;

        CommBufShift[0] = 0;
        CommBufShift[1] = localYSize;
        CommBufShift[2] = 2 * localYSize;
        CommBufShift[3] = 2 * localYSize + localXSize;
        CommBufShift[4] = 2 * (localXSize + localYSize);
        CommBufShift[5] = 2 * (localXSize + localYSize) + 1;
        CommBufShift[6] = 2 * (localXSize + localYSize) + 2;
        CommBufShift[7] = 2 * (localXSize + localYSize) + 3;

        CommBufSize[0] = localYSize;
        CommBufSize[1] = localYSize;
        CommBufSize[2] = localXSize; 
        CommBufSize[3] = localXSize;
        CommBufSize[4] = 1; 
        CommBufSize[5] = 1; 
        CommBufSize[6] = 1; 
        CommBufSize[7] = 1;
    }

    void ExchangeHalo(double* local) {
        PackBuf(local, sendBuf);
        vector<MPI_Request> req(16);
        int reqCnt = 0;

        for (int i = 0; i < 8; i++) {
            MPI_Isend(&sendBuf[CommBufShift[i]], CommBufSize[i], MPI_DOUBLE, Neighbor[i], i, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
            MPI_Irecv(&recvBuf[CommBufShift[i]], CommBufSize[i], MPI_DOUBLE, Neighbor[i], Tag[i], MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
        }

        MPI_Waitall(reqCnt, req.data(), MPI_STATUSES_IGNORE); 
        UnpackBuf(local, recvBuf); 
    }

    void PackBuf(const double* local, double* buf) {
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[0]; i++) {
            buf[CommBufShift[0] + i] = local[GetPos(0, i)];
        }
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[1]; i++) {
            buf[CommBufShift[1] + i] = local[GetPos(localXSize - 1, i)];
        }
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[2]; i++) {
            buf[CommBufShift[2] + i] = local[GetPos(i, 0)];
        }
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[3]; i++) {
            buf[CommBufShift[3] + i] = local[GetPos(i, localYSize - 1)];
        }

        buf[CommBufShift[4]] = local[GetPos(0, 0)];
        buf[CommBufShift[5]] = local[GetPos(0, localYSize - 1)];
        buf[CommBufShift[6]] = local[GetPos(localXSize - 1, 0)];
        buf[CommBufShift[7]] = local[GetPos(localXSize - 1, localYSize -1)];
    }

    void UnpackBuf(double* local, const double* buf) {
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[0]; i++) {
           local[GetSPos(0, i + 1)] = buf[CommBufShift[0] + i];
        }
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[1]; i++) {
            local[GetSPos(storedXSize - 1, i + 1)] = buf[CommBufShift[1] + i];
        }
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[2]; i++) {
            local[GetSPos(i + 1, 0)] = buf[CommBufShift[2] + i];
        }
        #pragma omp parallel for
        for (int i = 0; i < CommBufSize[3]; i++) {
            local[GetSPos(i + 1, storedYSize - 1)] = buf[CommBufShift[3] + i];
        }
        local[GetSPos(0, 0)] = buf[CommBufShift[4]];
        local[GetSPos(0, storedYSize - 1)] = buf[CommBufShift[5]];
        local[GetSPos(storedXSize - 1, 0)] = buf[CommBufShift[6]];
        local[GetSPos(storedXSize - 1, storedYSize - 1)] = buf[CommBufShift[7]];
    }

    int PlusMod(int val, int d, int mod) {
        int res = val + d + mod;
        while (res >= mod) {
            res -= mod;
        }
        return res;
    }

    int GetDest(int dx, int dy) {
        int x = PlusMod(mpiX, dx, mpiN);
        int y = PlusMod(mpiY, dy, mpiM);
        return GetMpiRank(x, y);
    }

    int GetMpiRank(int x, int y) {
        return mpiM * x + y;
    }

    int GetPos(int i, int j) {
        return GetSPos(i + 1, j + 1);
    }

    int GetSPos(int i, int j) {
        return storedYSize * i + j;
    }

    int GetX(int i) {
        return i + xStart;
    }

    int GetY(int j) {
        return j + yStart;
    }

    void InitialCondition(double* v1, double* v2) {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < storedXSize; i++) {
            for (int j = 0; j < storedYSize; j++) {
                int pos = GetSPos(i, j);
                double x = hx * GetX(i - 1);
                double y = hy * GetY(j - 1);
                v1[pos] = Phi1(x, y);
                v2[pos] = Phi2(x, y);
            }
        }
    }

    void InitV1x(double* v) {
        if (mpiX == 0) {
            #pragma omp parallel for
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(0, j)] = 0;
            }
        }
        if (mpiX == mpiN - 1) {
            #pragma omp parallel for
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(localXSize - 1, j)] = 0;
            }
        }
    }

    void ExchangeCyclicV1(double* local) {
        vector<MPI_Request> req(4);
        int reqCnt = 0;

        if (mpiY == mpiM - 1) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                cyclicSendBuf[i] = local[GetPos(i, localYSize - 1)];
            }
            int rightNeighbor = GetDest(0, 1);
            MPI_Isend(cyclicSendBuf, localXSize, MPI_DOUBLE, rightNeighbor, 1, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
            MPI_Irecv(cyclicRecvBuf, localXSize, MPI_DOUBLE, rightNeighbor, 0, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;

        }

        if (mpiY == 0) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                cyclicSendBuf[localXSize + i] = local[GetPos(i, 1)];
            }
            int leftNeighbor = GetDest(0, -1);
            MPI_Isend(&cyclicSendBuf[localXSize], localXSize, MPI_DOUBLE, leftNeighbor, 0, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
            MPI_Irecv(&cyclicRecvBuf[localXSize], localXSize, MPI_DOUBLE, leftNeighbor, 1, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
        }

        MPI_Waitall(reqCnt, req.data(), MPI_STATUSES_IGNORE); 

        if (mpiY == mpiM - 1) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                local[GetPos(i, localYSize)] = cyclicRecvBuf[i];
            }
        }

        if (mpiY == 0) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                local[GetPos(i, 0)] = cyclicRecvBuf[localXSize + i];
            }
        }
    }

    void ExchangeCyclicV2(double* local) {
        vector<MPI_Request> req(2);
        int reqCnt = 0;

        if (mpiY == 0) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                cyclicSendBuf[i] = local[GetPos(i, 1)];
            }
            int leftNeighbor = GetDest(0, -1);
            MPI_Isend(cyclicSendBuf, localXSize, MPI_DOUBLE, leftNeighbor, 0, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
        }

        if (mpiY == mpiM - 1) {
            int rightNeighbor = GetDest(0, 1);
            MPI_Irecv(cyclicRecvBuf, localXSize, MPI_DOUBLE, rightNeighbor, 0, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
        }

        MPI_Waitall(reqCnt, req.data(), MPI_STATUSES_IGNORE); 

        if (mpiY == mpiM - 1) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                local[GetPos(i, localYSize)] = cyclicRecvBuf[i];
            }
        }
    }

    void InitV1y(double* next) {
        ExchangeCyclicV1(next);
    }

    void InitV2x(double* v) {
        if (mpiX == 0) {
            #pragma omp parallel for
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(0, j)] = (4.0 * v[GetPos(1, j)] - v[GetPos(2, j)]) / 3.0;
            }
        }
        if (mpiX == mpiN - 1) {
            #pragma omp parallel for
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(localXSize - 1, j)] = (4.0 * v[GetPos(localXSize - 2, j)] - v[GetPos(localXSize - 3, j)]) / 3.0;
            }
        }
    }

    void InitV2y(double* v) {
        if (mpiY == 0) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                v[GetPos(i, 0)] = (4.0 * v[GetPos(i, 1)] - v[GetPos(i, 2)]) / 3.0;
            }
        }

        if (mpiY == mpiM - 1) {
            #pragma omp parallel for
            for (int i = 0; i < localXSize; i++) {
                v[GetPos(i, localYSize - 1)] = (4.0 * v[GetPos(i, localYSize - 2)] - v[GetPos(i, localYSize - 3)]) / 3.0;
            }
        }

        ExchangeCyclicV2(v);
    }

public:

    Solution(int argc, char** argv) {

        double time = MPI_Wtime();
        Init(argc, argv);

        InitialCondition(v1[0], v2[0]);

        for (int it = 0; it < gridT - 1; it++) {
            MPI_Barrier(MPI_COMM_WORLD);

            int cur = it % 2;
            int nxt = (it + 1) % 2;

            #pragma omp parallel for collapse(2)
            for (int i = 0; i < localXSize; i++) {
                for (int j = 0; j < localYSize; j++) {
                    int y = GetY(j);
                    int x = GetX(i);
                    {
                        bool edgeY = (y == 0);
                        bool edgeX = (x == 0 || x == gridN - 1);

                        // y == gridM - 1 is not an edge we have cyclic condition here
                        if (edgeX || edgeY) {
                            continue;
                        }

                        double x = GetX(i) * hx;
                        double y = GetY(j) * hy;
                        double t = ht * it;

                        Calc iterCalc(localYSize, hx, hy, ht, i, j, x, y, t);
                        iterCalc.CalcNextV1(v1[nxt], v1[cur], v2[cur]);
                        iterCalc.CalcNextV2(v2[nxt], v1[cur], v2[cur]);
                    }
                }
            }

            InitV1x(v1[nxt]);
            InitV1y(v1[nxt]);
            InitV2x(v2[nxt]);
            InitV2y(v2[nxt]);

            ExchangeHalo(v1[nxt]);
            ExchangeHalo(v2[nxt]);
           
            ExchangeCyclicV1(v1[nxt]);
            ExchangeCyclicV1(v2[nxt]);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        time = MPI_Wtime() - time;

        double delta1 = 0, delta2 = 0;
        for (int i = 0; i < localXSize; i++) {
            for (int j = 0; j < localYSize; j++) {
                int xi = GetX(i);
                int yj = GetY(j);
                double x = hx * xi;
                double y = hy * yj;
                double t = ht * (gridT - 1);
                int pos = GetPos(i, j);
                double val1 = u1(x,y,t);
                double val2 = u2(x,y,t);
                int it = (gridT - 1) % 2;
                double diff1 = abs(val1 - v1[it][pos]);
                double diff2 = abs(val2 - v2[it][pos]);
                delta1 += diff1 * diff1;
                delta2 += diff2 * diff2;
                if (printDebug)
                    cout << diff1 << " ";
            }
            if (printDebug) {
                cout << endl;
            }
        }
        if (printDebug) {
            cout << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double deltaMax1 = 0;
        MPI_Reduce(&delta1, &deltaMax1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        double deltaMax2 = 0;
        MPI_Reduce(&delta2, &deltaMax2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        printf("RES %d time %le delta1 %le delta2 %le\n", mpiRank, time, sqrt(hx * hy * delta1), sqrt( hx * hy * delta2));
        if (mpiRank == 0) {
            printf("RES Global wallTime %le ht %le delta1 %le delta2 %le sum %le\n", time, ht, sqrt(hx * hy * deltaMax1), sqrt(hx * hy * deltaMax2), sqrt( hx * hy * (deltaMax1 + deltaMax2)));
        }
    }

    // t = 0 Initial condition
    double Phi1(double x, double y) {
       return u1(x, y, 0);
    }

    // t = 0 Initial condition
    double Phi2(double x, double y) {
       return u2(x, y, 0);
    }

};


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Solution(argc, argv);
    MPI_Finalize();
}
