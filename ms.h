#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define inline __inline  // for VC
#define min __min  // for GCC
#define max __max

#define H 16
#define W 30
#define MINE 99
#define REPEAT 100
#define SAFE_DIST 0   // 开局点击格周围不放雷的距离(切比雪夫距离)，XP=0，WIN7=1
#define WHITE_HOLE 0  // 模拟每日挑战中间的不可点区域，2即4x4

// 这仨得是等差，且公差大于8，方便(un)markSharedTiles修改，也方便trim记录周围已知雷数
#define CLICK 10      // 可点格，之后会加上周围已知雷数，或标为公共区和非公共区
#define INDEP 20      // 格周围只有一个数字，不属于任意两数字公共区
#define SHARE 30      // 格周围至少有两个大于零数字

#define ISCLICK(x) ((x) >= CLICK && (x) < INDEP)
#define ISINDEP(x) ((x) >= INDEP && (x) < SHARE)
#define ISSHARE(x) ((x) >= SHARE && (x) < SHARE + (SHARE-INDEP))

#define show_error(msg) { printf("%s\n", msg); output(board); output(trim); return -1; }

// 要初始化必须用static
static int dir[][2] = {{-1, 0}, { 1,0}, {0,-1}, {0,1},  // ↑ ↓ ← →
                       {-1,-1}, {-1,1}, {1,-1}, {1,1}}; // ↖ ↗ ↙ ↘


// 定义时参数比声明这少(调用与声明一致)不报错，要注意
void output(int*);
void markSharedTiles(int*, int, int);
void unmarkSharedTiles(int*, int, int);
void checkSurroundTiles(int*, int, int, int*);

int click(int*, int*, int, int);
int mine_trim(int*, int*, int, int);
int logic_kernel(int*, int*, int);
int guess_kernel(int*, int*, int, double*, int*);


// 在头文件定义必须inline或static，否则链接时报错重定义
static inline int outofbound(int i, int j, int k)
{ return (i+dir[k][0]<0 || i+dir[k][0]>=H || j+dir[k][1]<0 || j+dir[k][1]>=W); }

static inline int sub2idx(int row, int col)
{ return row*W + col; }

static inline void idx2sub(int idx, int *row, int *col)
{ *row = idx / W;  *col = idx % W; }  // 快速取余A % M = A - A/M*M = idx - *row*W降低并行度?


#if H*W < MINE + (SAFE_DIST*2 + 1)*(SAFE_DIST*2 + 1) + WHITE_HOLE*WHITE_HOLE*4
#pragma message("TOO MANY MINES!")
MINE = 0;  // 弄个错误停止编译
#endif
