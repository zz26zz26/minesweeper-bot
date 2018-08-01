#include "ms.h"
#include <time.h>
#include <omp.h>

extern int dir[][2];


void output(int *board)
{
    int i, j, t, o = 0, m = 0;
    //system("cls");
    for (i = 0; i < H; i++) {
        for (j = 0; j < W; j++) {
            t = board[sub2idx(i, j)];
            o += (t < 0);
            m += (t == -9);
            switch (t) {
                case -99:  printf(" #"); break;  // 不可点击格
                case -10:  printf(" o"); break;  // 已点开的空白格
                case  -9:  printf(" $"); break;  // 已发现的雷(旗)
                case   0:  printf(" ."); break;  // 未点开的空白格
                case   9:  printf(" x"); break;  // 未发现的雷
                case  10:  printf(" ?"); break;  // 未知格(泛指 trim初始化为此)
                case  99:  printf("XX"); break;  // 死于此雷
                default:
                    if (t > CLICK)
                        printf("%d?", t % CLICK);
                    else if (t < -50)      // 在click过程中被标记为要点开
                        printf("-?");
                    else                   // 绝对值1-8直接输出；负数表示已点开
                        printf("%2d", t);  // 正数临时记录未点开格的周围已知雷数
            }
        }
        printf("\n");
    }
    printf("open = %d, mines = %d\n", o, m);
    system("pause");

    // trim导入matlab时，先把负号换成空格，再分别换? $ o为-1 -2 0
    // csv导入后board=getBoundedBoard(Untitled);
}

void print_progress(int s, clock_t t)
{
    int r;
    double v, dt;
    const double a = 0.250;                 // 新速度权重，越大越接近新值
    static double prev_v;                   // 打印是单线程
    static int prev_s;                      // static自动初始化0
    static clock_t prev_t;
    dt = (t - prev_t) / (double)CLOCKS_PER_SEC;
    if (dt > 1) {                           // 1s才更新一次各种数据
        v = (s - prev_s) / dt;              // 这段时间的速度，总平均也不很稳
        v = (prev_v > 0) ? (1 - a)*prev_v + a*v : v;
        r = (REPEAT - s) / v;
        prev_s = s;
        prev_t = t;
        prev_v = v;
    }
    else {
        v = prev_v;
        r = (REPEAT - s) / v;               // 别的时候也会用之前速度重算剩余
    }

    printf("%6.2f%%", s * 100.0 / REPEAT);  // 覆写的不足6字符就换行会留下残余
    if (t > CLOCKS_PER_SEC) {
        //printf("%10.1f/s", v*H*W);
        if (r > 5400)     printf("%4dh  \b\b\b\b\b\b\b", (r + 3599) / 3600);  // 2.1要显示3；参考复制5s为单位
        else if (r > 90)  printf("%4dmin\b\b\b\b\b\b\b", (r >= 600) ? (r + 299) / 300 * 5 : (r + 59) / 60);
        else              printf("%4ds  \b\b\b\b\b\b\b", (r >= 10) ? (r + 4) / 5 * 5 : r + 1);
        //printf("\b\b\b\b\b\b\b\b\b\b\b\b");
    }
    printf("\b\b\b\b\b\b\b");
}

// 0周围的格子加入待点开队列，由于0周围不会有雷，周围所有没点开的入队
// trim中的0代表周围没有未知的雷，亦可把所有没点开的入队
// 另外利用trim把要点开但还没点开的-50标记为负数不重复入队，同时可还原记录的已知雷数
#define push_tiles_around_0 {                                   \
    idx2sub(r, &i, &j);                                         \
    for (k = 0; k < 8; k++) {                                   \
        if (outofbound(i, j, k))  continue;                     \
        s = sub2idx(i + dir[k][0], j + dir[k][1]);              \
        if (trim[s] >= 0) { trim[s] -= 50; num[tail++] = s; }   \
    }                                                           \
} // 写成函数需参数(board/r/num/tail)，内部变量(i/j/k/s)

int click(int *board, int *trim, int row, int col)
{
    int i, j, k, r, s, area, head, tail, num[H*W];

    area = 0;             // 按8连通区域扩展面积
    head = tail = 0;      // 没重复，不必循环队列
    r = num[tail++] = sub2idx(row, col);

    if (board[r] == 9) {  // 见雷死，死个明白
        board[r] = trim[r] = 99;
        return -1;
    }

    while (head != tail) {         // BFS, flood-fill
        r = num[head++];           // pop，正宗队列，减少重复访问

        if (board[r] == 0) {
            board[r] = -10;        // mark visited
            area++;
            push_tiles_around_0;
            trim[r] = board[r];    // push不管中心格；覆写不必+100
        }
        else if (board[r] > 0) {   // 正数只有数字(0-8)会入队
            board[r] = -board[r];  // mark visited
            area++;

            // 调用click时trim中的INDEP、H_MINE等已还原，最大10最小-10，push中减去后应小于原最小值
            trim[r] = board[r] + (trim[r] - CLICK);
            if (trim[r] > -90 && trim[r] < -10) {   // board中除了障碍最小-8，此条件筛出push中标记的待点开格
                trim[r] += 50;                      // 去掉"要点开"标记(注意别搞没标记就点开的，如参数指定)
            }
            if (trim[r] == 0) {
                trim[r] = -10;
                push_tiles_around_0;
            }
        }
        // else即(board[r] < 0)表示访问过，不增加area
        //output(trim);
    }
    return area;
}

// 异或交换不太用寄存器(寄存器传的参是地址)，且相等的数异或得0
void swap(int *a, int *b)
{
    int tem = *a;
    *a = *b;
    *b = tem;
}

// 点击指定位置开局，并返回展开区域面积
int open_new(int *board, int *trim, int row, int col)
{
    int i, j, k, m, r, s, len, num[H*W];

    len = H*W;
    memset(board, 0, sizeof(int)*len);               // memset单位是char，不是int
    memcpy(trim, board, sizeof(int)*len);            // trim里未知格统一设为可点，不含其他信息
    for (k = 0; k < len; k++)  trim[k] = CLICK;      // 初始必有board[k]>=0，不必if

    for (k = 0; k < len; k++)  num[k] = k;           // 生成不重复随机数模拟放雷
    swap(&num[r = sub2idx(row, col)], &num[--len]);  // 排除点的第一格，那里不能放雷

#if SAFE_DIST + WHITE_HOLE > 0                       // win7后还会把点第一格周围8格也排除在可放雷位置外
    int hash[H*W];                                   // 这个开内存不很耗时
    for (k = 0; k < H*W; k++)  hash[k] = k;          // len减过1了不能再用k<len了啊
    swap(&hash[num[r]], &hash[num[len]]);            // 中心格对应前面先换，不是换hash的[r]和[len]啊

    // hash[i]记录最初值是i的格现在的下标，跟着数字一起换
    // 每次把一个该换的换到位，从num末尾往前排(中心格第一个排)，有时会交换俩该换的 但最后会把可用的换到前面
#if SAFE_DIST > 0
    for (i = row - SAFE_DIST; i <= row + SAFE_DIST; i++) {
        if (i < 0 || i >= H)  continue;
        for (j = col - SAFE_DIST; j <= col + SAFE_DIST; j++) {
            if (j < 0 || j >= W)  continue;
            if (i == row && j == col)  continue;
            s = sub2idx(i, j);
            swap(&num[hash[s]], &num[--len]);
            swap(&hash[num[hash[s]]], &hash[num[len]]);
        }
    }
#endif

#if WHITE_HOLE > 0
    for (i = H / 2 - WHITE_HOLE; i < H / 2 + WHITE_HOLE; i++) {
        if (i < 0 || i >= H)  continue;
        for (j = W / 2 - WHITE_HOLE; j < W / 2 + WHITE_HOLE; j++) {
            if (j < 0 || j >= W)  continue;
            s = sub2idx(i, j);
            board[s] = trim[s] = -99;
            swap(&num[hash[s]], &num[--len]);
            swap(&hash[num[hash[s]]], &hash[num[len]]);
        }
    }
#endif
#endif

    for (m = 0; m < MINE; m++) {
        r = rand() % len;
        for (k = 0; k < 8; k++) {
            idx2sub(num[r], &i, &j);
            if (outofbound(i, j, k))  continue;
            s = sub2idx(i + dir[k][0], j + dir[k][1]);
            if (board[s] >= 0 && board[s] < 9)  // 不能改动附近的雷或空洞的格
                board[s]++;                     // 雷周围格子里的数字+1
        }
        board[num[r]] = 9;
        num[r] = num[--len];  // 用过直接排除到可放雷位置外，类似抽签不放回仍为等概率
    }

    return click(board, trim, row, col);
}

// loop: trim -> guess -> click
// board用于存当前进度和雷位置(已知格<0，未点开>=0)
// trim存由board推出的每个数字旁已知雷数，累加到未知可点格和已知数字(见mine_trim)的值上，因此也有>0的格
// 数字变化以雷为中心，推理时找到一个雷就把他周围数字值都-1，减到0就点开周围可点格，省去最后全盘trim再点
int play_game(int *board, int row, int col)
{
    int area, steps, mines, num[H*W], trim[H*W];
    double prob[H*W];

    mines = 0;
    steps = 1;
    area = open_new(board, trim, row, col);
    //output(board);
    //return area;  // open_new不初始化trim仅快5%

    while (1) {
        // infer on 'board', trim, click
        mines = logic_kernel(board, trim, mines);
        if (mines == MINE)  break;  // WIN!!! 但存在3x3的雷区时仅用logic会死循环
        if (mines < 0)  return 0;   // DIE
        //printf("logic finish:\n"); output(trim);

        // guess on 'trim' then find min prob to click
        mines = guess_kernel(board, trim, mines, prob, num);
        if (mines == MINE)  break;  // WIN
        if (mines < 0)  return 0;   // DIE
    }

    return steps;
}


int main()
{
    clock_t a, b;
    double sum = 0;
    int i, j, n, k = 0;
    int board[H*W], count[H*W] = { 0 };

    a = clock();

    #pragma omp parallel private(i, j, n, board) //num_threads(1)
    {
        srand(time(NULL) + omp_get_thread_num());  // 多线程的随机数得注意一下
        srand(rand() - rand() + rand());           // 第一个随机数接近seed；参数无符号，负数没事

        // 遇到omp parallel，就会创造多个线程跑代码(循环也有多个在跑)
        // 但见omp for，又会把循环拆块分给这些线程跑；但用parallel for嵌套一般会慢
        // 并且omp atomic只支持标量算数运算，几乎不拖慢速度；critical放内层循环就比较明显了
        // 另外omp single/for/sections同属工作分配构造；master/atomic等属同步构造；前者不可嵌套
        // 最后omp parallel/for/sections/single末尾含隐式barrier(可在omp for之类后加nowait去掉)
        // PS：如果board为全局还可在声明后用omp threadprivate指定为线程全局，但不能用动态线程
        // PPS：C++才能在omp for里用for(int i...)；另发现改x64能快15%
        #pragma omp for schedule(dynamic)
        for (n = 0; n < REPEAT; n++) {
            for (i = 0; i < H; i++) {
                for (j = 0; j < W; j++) {
                    #pragma omp atomic
                    count[sub2idx(i, j)] += play_game(board, i, j);
                    //printf("tid=%d, (%2d,%2d,%2d))\n", omp_get_thread_num(), n, i, j);
                }   //printf("tid=%d, %d\n", omp_get_thread_num(), rand());
            }

            #pragma omp critical  // 输出进度，初级约慢10%，加的快不用atomic竟然还会不一致
            {
                #pragma omp atomic
                k++;
                print_progress(k, clock() - a);
            }
        }

        #pragma omp for reduction(+: sum)
        for (n = 0; n < H*W; n++) {  // k是shared
            sum += count[n];
        }
    }

    b = clock();

    FILE *out = fopen("stat.csv", "w");
    if (out == NULL)  out = stdout;
    for (i = 0; i < H; i++) {
        for (j = 0; j < W; j++)
            fprintf(out, "%d,", count[sub2idx(i, j)]);
        fprintf(out, "\n");
    }
    fclose(out);

    printf("wins = %.3f%%\n", sum * 100.0 / (H*W*REPEAT));  // 下面不转可能溢出int；别cls 错误信息没了
    printf("time = %dms (%.2f game/s)  ", b - a, ((double)H*W*REPEAT*CLOCKS_PER_SEC) / (b - a));
    system("pause");
    return 0;
}
