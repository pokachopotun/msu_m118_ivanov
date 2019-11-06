#include "mpi.h"
#include <iostream>
#include <vector>

using namespace std;

void PackBuf(const double* local, int localXSize, int localYsize, double* buf) {
    int storedXSize = localXSize + 2;
    int storedYSize = localYSize + 2;
    int n = localXSize * localYSize;
    int N = storedXSize * storedYSize;

    int shift = 0;
    for (int i = 0; i < localYSize; i++) {
        buf[shift + i] = local[storedYSize + 1 + i];
    }
    shift += localYSize;
    for (int i = 0; i < localYSize; i++) {
        buf[shift + i] = local[storedYSize * (storedXSize - 2) + 1 + i];
    }
    shift += localYSize;
    for (int i = 0; i < localXSize; i++) {
         buf[shift + i] = local[storedYSize * (i + 1) + 1];
    }
    shift += localXSize;
    for (int i = 0; i < localXSize; i++) {
        buf[shift + i] = local[storedYSize * (i + 2) - 2];
    }
    shift += localXSize;
    buf[shift++] = local[storedYSize + 1];
    buf[shift++] = local[2 * storedYSize - 2];
    buf[shift++] = local[N - 2 * storedYSize + 1];
    buf[shift++] = local[N - storedYSize - 2];
}

void UnpackBuf(double* local, int localXSize, int localYsize, const double* buf) {
    int storedXSize = localXSize + 2;
    int storedYSize = localYSize + 2;
    int n = localXSize * localYSize;
    int N = storedXSize * storedYSize;

    int shift = 0;
    for (int i = 0; i < localYSize; i++) {
       local[storedYSize + 1 + i] = buf[shift + i];
    }
    shift += localYSize;
    for (int i = 0; i < localYSize; i++) {
        local[storedYSize * (storedXSize - 2) + 1 + i] = buf[shift + i];
    }
    shift += localYSize;
    for (int i = 0; i < localXSize; i++) {
        local[storedYSize * (i + 1) + 1] = buf[shift + i];
    }
    shift += localXSize;
    for (int i = 0; i < localXSize; i++) {
        local[storedYSize * (i + 2) - 2] = buf[shift + i];
    }
    shift += localXSize;
    local[storedYSize + 1] = buf[shift++];
    local[2 * storedYSize - 2] = buf[shift++];
    local[N - 2 * storedYSize + 1] = buf[shift++];
    local[N - storedYSize - 2] = buf[shift++];
}

bool CheckRank(int mpiNodes, int mpiRank) {
    return 0 <= mpiRank && mpiRank < mpiNodes;
}

void DoSend(double* sendBuf, int shift, int cnt, int dest, int mpiNodes, int vector<MPI_Request>& req, int& reqCnt) {
    if (CheckRank(mpiNodes, dest)) {
        MPI_Isend(sendBuf[shift], cnt, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &req[reqCnt]);
        reqCnt++;
    }
}

void DoRecv(double* recvBuf, int shift, int cnt, int dest, int mpiNodes, int vector<MPI_Request>& req, int& reqCnt) {
    if (CheckRank(mpiNodes, dest)) {
        MPI_Irecv(recvBuf[shift], cnt, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &req[reqCnt]);
        reqCnt++;
    }
}

void ExchangeHalo(double* local, double* sendBuf, vector<double>& recvBuf, int localXSize, int localYSize, int mpiN, int mpiM, int mpiRank) {
    PackBuf(local, localXSize, localYSize, sendBuf);
    vector<MPI_Request> req(32);
    reqCnt = 0;

    int mpiNodes = mpiN * mpiM;
    int shift = 0;

    //send top
    DoSend(sendBuf, 0, localYSize, mpiRank - mpiN, mpiNodes, req, reqCnt);
    
    //send bottom
    DoSend(sendBuf, localYSize, localYSize, mpiRank + mpiN, mpiNodes, req, reqCnt);

    //send left
    DoSend(sendBuf, 2 * localYSize, localXSize, mpiRank - 1, mpiNodes, req, reqCnt);

    //send right
    DoSend(sendBuf, 2 * localYSize + localXSize, localYSize, mpiRank + 1, mpiNodes, req, reqCnt);

    //send top left
    DoSend(sendBuf, 2 * (localXSize + localYSize), 1, mpiRank - mpiN - 1, mpiNodes, req, reqCnt);

    //send top right
    DoSend(sendBuf, 2 * (localXSize + localYSize) + 1, 1, mpiRank - mpiN + 1, mpiNodes, req, reqCnt);

    //send bottom left
    DoSend(sendBuf, 2 * (localXSize + localYSize) + 2, 1, mpiRank + mpiN - 1, mpiNodes, req, reqCnt);

    //send bottom right
    DoSend(sendBuf, 2 * (localXSize + localYSize) + 3, 1, mpiRank + mpiN + 1, mpiNodes, req, reqCnt);

    //recv top
    DoRecv(recvBuf, 0, localYSize, mpiRank - mpiN, mpiNodes, req, reqCnt);
    
    //recv bottom
    DoRecv(recvBuf, localYSize, localYSize, mpiRank + mpiN, mpiNodes, req, reqCnt);

    //recv left
    DoRecv(recvBuf, 2 * localYSize, localXSize, mpiRank - 1, mpiNodes, req, reqCnt);

    //recv right
    DoRecv(recvBuf, 2 * localYSize + localXSize, localYSize, mpiRank + 1, mpiNodes, req, reqCnt);

    //recv top left
    DoRecv(recvBuf, 2 * (localXSize + localYSize), 1, mpiRank - mpiN - 1, mpiNodes, req, reqCnt);

    //recv top right
    DoRecv(recvBuf, 2 * (localXSize + localYSize) + 1, 1, mpiRank - mpiN + 1, mpiNodes, req, reqCnt);

    //recv bottom left
    DoRecv(recvBuf, 2 * (localXSize + localYSize) + 2, 1, mpiRank + mpiN - 1, mpiNodes, req, reqCnt);

    //recv bottom right
    DoRecv(recvBuf, 2 * (localXSize + localYSize) + 3, 1, mpiRank + mpiN + 1, mpiNodes, req, reqCnt);

    MPI_Waitall(reqCnt, req.data(), MPI_STATUSES_IGNORE); 
    UnpackBuf(local, localXSize, localYSize, recvBuf); 
}

class Calc {
private:
    int tl;
    int t;
    int tr;

    int l;
    int m;
    int r;

    int bl;
    int b;
    int br;

    double hx;
    double hy;

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
    Calc(int localYSize, double hx, double hy, int i, int j) : hx(hx), hy(hy) {
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
       
        next1[m] = v1[m] + laplas1 + g1;
        next2[m] = v2[m] + laplas2 + g2;
    }
};

int main(int argc, char** argv[]) {
    MPI_Init(&argc, &argv);

    double xmax, ymax, tmax;

    int mpiN;
    int mpiM;
    int gridN;
    int gridM;
    int TN;

    int mpiNodes, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (mpiNodes != mpiN * mpiM) {
        cerr << "mpiNodes != mpiN * mpiM" << endl;
        return 1;
    }

    int mpiX = mpiRank / mpiN;
    int mpiY = mpiRank % mpiN;

    bool isXMax = curMpiX == mpiN - 1;  
    bool isYMax = curMpiY == mpiM - 1;

    int localXSize = isXMax ? gridN % mpiN : gridN / mpiN;
    int localYSize = isYMax ? gridM % mpiM : gridM / mpiM;
    int storedXSize = localXSize + 2;
    int storedYSize = localYSize + 2;

    vector<double*> v1(TN);
    vector<double*> v2(TN);

    int commBufSize = 2 * (localXSize + localYSize) + 4;

    double* sendBuf;
    double* recvBuf;

    //TODO: 1) zero and initial conditions. Fill matrix;
    //TODO: 2) init all the arrays
    //TODO: 3) calc answer diff

    // do iter
    for (int it = 0; it < TN - 1; it++) {
        
        ExchangeHalo(v1[it], sendBuf, recvBuf, localXSize, localYSize, mpiN, mpiM, mpiRank);
        ExchangeHalo(v2[it], sendBuf, recvBuf, localXSize, localYSize, mpiN, mpiM, mpiRank);

        for (int i = 0; i < localXSize; i++) {
            for (int j = 1; j < localYSize; j++) {
                auto iterCalc = Calc(localYSize, hx, hy, i, j);
                iterCalc.CalcNextT(v1[it + 1], v2[it + 1], v1[it], v2[it]);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
}
