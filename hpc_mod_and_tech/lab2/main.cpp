#include "mpi.h"
#include <iostream>
#include <vector>

using namespace std;


class Calc {
private:
    int tl, t, tr;
    int l, m, r;
    int bl, b, br;
    double hx, hy, ht;
    int i, j;

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
    Calc(int localYSize, double hx, double hy, double ht,  int i, int j) : hx(hx), hy(hy), ht(ht), i(i), j(j) {
        int storedYSize = localYSize + 2; 

        tl = storedYSize * (i - 1) + j;
        t = storedYSize * (i - 1) + 1 + j;
        tr = storedYSize * (i - 1) + 1 + j + 1;

        l = storedYSize * (i + 1) + j;
        m = storedYSize * (i + 1) + 1 + j;
        r = storedYSize * (i + 1) + 1 + j + 1;

        bl = storedYSize * (i + 2) + j;
        b = storedYSize * (i + 2) + 1 + j;
        br = storedYSize * (i + 2) + 1 + j + 1;
    }

    void CalcNextT(double* next1, double* next2, double* v1, double* v2) {
        double laplas1 = Laplas(v1);      
        double laplas2 = Laplas(v2);
        double g1 = Grad1(v1, v2);
        double g2 = Grad2(v1, v2);
       
        next1[m] = ht * (v1[m] + laplas1) + g1;
        next2[m] = ht * (v2[m] + laplas2) + g2;
    }
};

// t = 0 Initial condition
double Phi1(double x, double y) {
   return x + y;
}

// t = 0 Initial condition
double Phi2(double x, double y) {
   return x + y;
}

// var 10
// v1(x = 0, x = xmax) = 0
// v1(y) periodic
// v2(x) 2nd
// v2(y) 2nd



class Solution {
private:
    // mpi grid size
    int mpiN, mpiM;
    // solution grid size 
    int gridN, gridM, gridT;
    // constraints
    double xmax, ymax, tmax;
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
    int CommBufShift[8];
    int CommBufSize[8];

    void Init(int argv, char** argv) {
        // mpi grid size
        mpiN = stoi(argv[1]);
        mpiM = stoi(argv[2]);

        // solution grid size 
        gridN = stoi(argv[3]);
        gridM = stor(argv[4]);
        gridT = stoi(argv[5]);

        // constraints
        xmax = atof(argv[6]);
        ymax = atof(argv[7]);
        tmax = arof(argv[8]);

        // grid step
        hx = xmax / gridN;
        hy = ymax / gridM;
        ht = tmax / gridT; 

        MPI_Comm_size(MPI_COMM_WORLD, &mpiNodes);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        if (mpiNodes != mpiN * mpiM) {
            cerr << "mpiNodes != mpiN * mpiM" << endl;
            return;
        }

        mpiX = mpiRank / mpiN;
        mpiY = mpiRank % mpiN;

        bool isXMax = curMpiX == mpiN - 1;  
        bool isYMax = curMpiY == mpiM - 1;

        localXSize = isXMax ? gridN % mpiN : gridN / mpiN;
        if (localXSize == 0) {
            localXSize = gridN / mpiN;
        }
        localYSize = isYMax ? gridM % mpiM : gridM / mpiM;
        if (localYSize == 0) {
            localYSize = gridM / mpiM;
        }
        storedXSize = localXSize + 2;
        storedYSize = localYSize + 2;

        sliceSize = storedXSize * storedYSize;
        commBufSize = 2 * (localXSize + localYSize) + 4;

        xStart = mpiX * localXSize;
        yStart = mpiY * localYsize;

        v1.reserve(gridT);
        v2.reserve(gridT);
        
        for (int i = 0; i < gridT; i++) {
            v1.push_back(new double[sliceSize]);
            v2.push_back(new double[sliceSize]);
        }

        sendBuf = new double[commBufSize];
        recvBuf = new double[commBufSize];

        Neighbor = {
            GetDest(-1, 0),
            GetDest(1, 0),
            GetDest(0, -1),
            GetDest(0, 1),
            GetDest(-1, -1),
            GetDest(-1, 1),
            GetDest(1, -1),
            GetDest(1, 1)
        };

        CommBufShift = {
            0,
            localYSize,
            2 * localYSize,
            2 * localYSize + localXSize,
            2 * (localXSize + localYSize),
            2 * (localXSize + localYSize) + 1,
            2 * (localXSize + localYSize) + 2,
            2 * (localXSize + localYSize) + 3
        };

        CommBufSize = {
            localYSize, localYSize,
            localXSize, localXSize,
            1, 1, 1, 1
        };
    }

    void ExchangeHalo(double* local) {
        PackBuf(local, sendBuf);
        vector<MPI_Request> req(32);
        reqCnt = 0;

        for (int i = 0; i < 8; i++) {
            MPI_Isend(sendBuf[CommBufShift[i]], CommBufSize[i], MPI_DOUBLE, Neighbor[i], 0, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
            MPI_Irecv(recvBuf[CommBufShift[i]], CommBufSize[i], MPI_DOUBLE, Neighbor[i], 0, MPI_COMM_WORLD, &req[reqCnt]);
            reqCnt++;
        }

        MPI_Waitall(reqCnt, req.data(), MPI_STATUSES_IGNORE); 
        UnpackBuf(local, recvBuf); 
    }

    void PackBuf(const double* local, double* buf) {
        int n = localXSize * localYSize;
        int N = storedXSize * storedYSize;

        for (int i = 0; i < CommBufSize[0]; i++) {
            buf[CommBufShift[0] + i] = local[storedYSize + 1 + i];
        }
        for (int i = 0; i < CommBufSize[1]; i++) {
            buf[CommBufShift[1] + i] = local[storedYSize * (storedXSize - 2) + 1 + i];
        }
        for (int i = 0; i < CommBufSize[2]; i++) {
            buf[CommBufShift[2] + i] = local[storedYSize * (i + 1) + 1];
        }
        for (int i = 0; i < CommBufSize[3]; i++) {
            buf[CommBufShift[3] + i] = local[storedYSize * (i + 2) - 2];
        }

        buf[CommBufShift[4]] = local[storedYSize + 1];
        buf[CommBufShift[5]] = local[2 * storedYSize - 2];
        buf[CommBufShift[6]] = local[N - 2 * storedYSize + 1];
        buf[CommBufShift[7]] = local[N - storedYSize - 2];
    }

    void UnpackBuf(double* local, const double* buf) {
        int n = localXSize * localYSize;
        int N = storedXSize * storedYSize;

        for (int i = 0; i < CommBufSize[0]; i++) {
           local[storedYSize + 1 + i] = buf[CommBufShift[0]];
        }
        for (int i = 0; i < CommBufSize[1]; i++) {
            local[storedYSize * (storedXSize - 2) + 1 + i] = buf[CommBufShift[1]];
        }
        for (int i = 0; i < CommBufSize[2]; i++) {
            local[storedYSize * (i + 1) + 1] = buf[CommBufShift[2]];
        }
        for (int i = 0; i < CommBufSize[3]; i++) {
            local[storedYSize * (i + 2) - 2] = buf[CommBufShift[3]];
        }
        local[storedYSize + 1] = buf[CommBufShift[4]];
        local[2 * storedYSize - 2] = buf[CommBufShift[5]];
        local[N - 2 * storedYSize + 1] = buf[CommBufShift[6]];
        local[N - storedYSize - 2] = buf[CommBufShift[7]];
    }

    bool CheckRank(int mpiNodes, int mpiRank) {
        return 0 <= mpiRank && mpiRank < mpiNodes;
    }

    int PlusMod(int val, int d, int mod) {
        int res = val + d + mod;
        while (res >= mod) {
            res -= mod;
        }
        return res;
    }

    int GetDest(int dx, int dy) {
        int x = PlusMod(mpiX, dx);
        int y = PlusMod(mpiY, dy);
        return GetMpiRank(x, y);
    }

    int MpiRank(int x, int y) {
        return mpiM * x + y;
    }

    int GetPos(int i, int j) {
        return storedYSize * (i + 1) + 1 + j;
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
        for (int i = 0; i < localXSize; i++) {
            for (int j = 0; j < localYSize; j++) {
                int pos = GetPos(i, j);
                double x = hx * GetX(i);
                double y = hy * GetY(j);
                v1[pos] = Phi1(x, y);
                v2[pos] = Phi2(x, y);
            }
        }
    }

    void InitV1x(double* v) {
        if (mpiX == 0) {
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(0, j)] = 0;
            }
        }
        if (mpiX == mpiN - 1) {
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(localXSize - 1, j)] = 0;
            }
        }
    }

    void InitV1y(double* next, double* v) {
        if (mpiY == 0) {
            for (int i = 0; i < localXSize; i++) {
                next[GetPos(i, 0)] = v[GetSPos(i, 0)];
            }
        }
        if (mpiY == mpiM - 1) {
            for (int i = 0; i < localXSize; i++) {
                next[GetPos(i, localYSize - 1)] = v[GetSPos(i, localYSize - 1)];
            }
        }
    }

    void InitV2x(double* v) {
        if (mpiX == 0) {
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(0, j)] = (v[GetPos(1, j)] - v[GetPos(2, j)]) / 3.0;
            }
        }
        if (mpiX == mpiM - 1) {
            for (int j = 0; j < localYSize; j++) {
                v[GetPos(localXSize - 1, j)] = (4.0 * v[GetPos(localXSize - 2, j)] - v[GetPos(localXSize - 3, j)]) / 3.0;
            }
        }
    }

    void InitV2y(double* v) {
        if (mpiY == 0) {
            for (int i = 0; i < localYSize; i++) {
                v[GetPos(i, 0)] = (v[GetPos(i, 1)] - v[GetPos(i, 2)]) / 3.0;
            }
        }
        if (mpiY == mpiM - 1) {
            for (int i = 0; i < localYSize; i++) {
                v[GetPos(i, localXSize - 1)] = (4.0 * v[GetPos(i, localYSize - 2)] - v[GetPos(i, localYSize - 3)]) / 3.0;
            }
        }
    }

public:

    Solution(int argc, char** argv) {
        Init(argc, argv);

        InitialCondition(v1[0], v2[0]);

        for (int it = 0; it < gridT - 1; it++) {
            MPI_Barrier(MPI_COMM_WORLD);
            
            ExchangeHalo(v1[it]);
            ExchangeHalo(v2[it]);

            for (int i = 0; i < localXSize; i++) {
                for (int j = 0; j < localYSize; j++) {
                    if (GetX(i) == 0 || GetX(i) == gridN - 1 ) {
                        continue;
                    }
                    if (GetY(j) == 0 || GetY(j) == gridM - 1) {
                        continue;
                    }
                    auto iterCalc = Calc(localYSize, hx, hy, ht, i, j);
                    iterCalc.CalcNextT(v1[it + 1], v2[it + 1], v1[it], v2[it]);
                }
            }

            InitV1x(v1[it]);
            InitV1y(v1[it + 1], v1[it]);
            InitV2x(v2[it + 1]);
            InitV2y(v2[it + 1]);
        }
    }
};    

int main(int argc, char** argv[]) {
    MPI_Init(&argc, &argv);
    
    //TODO: 3) calc answer diff

    Solution(argc, argv);

    MPI_Finalize();
}
