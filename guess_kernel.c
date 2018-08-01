#include "ms.h"
// fp:fast: https://msdn.microsoft.com/zh-cn/library/e7s85ffb.aspx

#define H_MINE 40   // 假设是雷，要与board原有的不同
#define H_SAFE 50   // 其余是假设是安全，不能留可点
#define I_PROB 2.0  // 没让猜的格子的概率，应大于1，且float可精确表示(如2的整数次幂)

typedef struct {
    int res_capacity;
    int res_length;
    int *res_hash;  // 合并case时的辅助数组

    int *num_mine;
    double *freq_board;
    double *case_count;
} res_info;  // 排列结果及其辅助数组

extern int dir[][2];  // 说明在文件外定义

// 输入的board为trim之后
// used - 在数字旁的可点格，用排列组合搞出每格概率
// rest - 不在数字旁的可点格，概率都一样
int guess_kernel(int *board, int *trim, int found_mine, double *prob, int *min_loc)
{
    int *r_idx, *c_idx, max_depth, max_mine, retVal;
    int num_tiles, rest_tile;
    int used_mine, rest_mine;
    double used_case, rest_case, sum_case;
    double *used_freq, rest_freq;
    double *freq_unshared, *inc_freq_board;
    res_info used;  // 在数字旁的可点格的各种概率

    int cnk(int, int);
    double Cnk(int, int);
    int try_mine_recur(int*, int*, int*, int, int, int, int, double,
                       double*, double*, res_info*);

    max_mine = MINE - found_mine;
    if (max_mine < 0)  return -1;  // ==0应该赢了，进来可以把被雷包围的可点格找出

    num_tiles = H * W;
    used.res_length = 1;
    used.res_capacity = 2;  // 初值至少2，增加时不必判断
    r_idx = (int*)malloc(sizeof(int)*num_tiles);  // 多点就多点了
    c_idx = (int*)malloc(sizeof(int)*num_tiles);
    used.res_hash   = (int*)malloc(sizeof(int)*(max_mine + 1));  // 要+1啊傻
    used.num_mine   = (int*)malloc(sizeof(int)*used.res_capacity);
    used.case_count = (double*)malloc(sizeof(double)*used.res_capacity);
    used.freq_board = (double*)malloc(sizeof(double)*num_tiles*used.res_capacity);
    freq_unshared   = (double*)malloc(sizeof(double)*num_tiles);
    inc_freq_board  = (double*)malloc(sizeof(double)*num_tiles);

    memset(used.res_hash,  -1, sizeof(int)*(max_mine + 1));
    memset(used.num_mine,   0, sizeof(int)*used.res_capacity);
    memset(used.case_count, 0, sizeof(double)*used.res_capacity);
    memset(used.freq_board, 0, sizeof(double)*num_tiles*used.res_capacity);
    memset(freq_unshared,   0, sizeof(double)*num_tiles);


    // max_depth还是0可能因为开局角落是3/一些格被已知雷包围，trim后没数字
    max_depth = 0;
    for (int k = 0; k < H*W; k++) {
        if (trim[k] > -9 && trim[k] < 0) {
            idx2sub(k, &r_idx[max_depth], &c_idx[max_depth]);
            markSharedTiles(trim, r_idx[max_depth], c_idx[max_depth]);
            max_depth++;
        }
    }

    // trim已trim过，且已知格为负数，在里面做了转换
    // 且不必再用realmin表示算过的格安全，数字旁可点肯定都是，0即可认为安全
    retVal = try_mine_recur(trim, r_idx, c_idx, max_mine, max_depth,
                            0, 0, 1.0, freq_unshared, inc_freq_board, &used);

    if (retVal == 0) {
        sum_case = 0;
        rest_tile = 0;
        memset(prob, 0, sizeof(double)*num_tiles);  // 后面是+=，得初始化啊啊啊
        for (int k = 0; k < H*W; k++) {  // 没被markShared的远离数字的可点格
            rest_tile += ISCLICK(trim[k]);
        }

        for (int s = 0; s < used.res_length; s++) {
            used_mine = used.num_mine[s];
            rest_mine = MINE - found_mine - used_mine;
            if (used_mine > max_mine || rest_mine > rest_tile)  continue;

            // 用double才够存组合数，注意函数名C大写
            rest_case = Cnk(rest_tile,   rest_mine);    // 开局角落1，C(476,98)>1e103
            rest_freq = Cnk(rest_tile-1, rest_mine-1);  // 每个数字在C(n,k)中出现C(n-1,k-1)次
            used_case = used.case_count[s];
            used_freq = &(used.freq_board[H*W*s]);
            for (int k = 0; k < H*W; k++) {  // 在数字旁和不在数字旁的未知格情况数交叉相乘
                if (ISCLICK(trim[k]))   // 不在数字旁的未知格(可点格)
                    prob[k] += rest_freq * used_case;
                else if (trim[k] >= 0)  // 数字旁未知格 (trim再mark后未知只剩可点/数字旁独立/公共)
                    prob[k] += used_freq[k] * rest_case;  // 数字旁可点肯定都是算过的
            }
            sum_case += used_case * rest_case;
        }

        for (int k = 0; k < H*W; k++) {  // 算总概率
            if (trim[k] >= 0)
                prob[k] /= sum_case;
            else
                prob[k] = I_PROB;  // 已知格(雷/数字)要设为大数字
        }
    }

    // 恢复现场让后面发现确定的雷时点开周围不出错
    // 注意上面会用到公共区等信息来区分数字周围的格，算完概率再恢复
    for (int k = 0; k < max_depth; k++) {
        unmarkSharedTiles(trim, r_idx[k], c_idx[k]);
    }

    // free(NULL)不会出错
    free(r_idx);
    free(c_idx);
    free(used.res_hash);
    free(used.num_mine);
    free(used.case_count);
    free(used.freq_board);
    free(freq_unshared);
    free(inc_freq_board);

    if (retVal != 0)  show_error("error in guess.");  // 排列出错就算输咯

    // 以下先对概率做点手脚，然后找到最小概率格点开
    int i, j, k, kk, clik, mins;
    double minprob, temprob[H*W] = { 0 };
    // 做手脚，多点边缘(不必考虑开局，因为open_new已经点开第一格了)
    for (k = 0; k < H*W; k++) {
        if (prob[k] >= 1.0)  continue;  // 一定是雷的别动
        temprob[k] = 1.0 - prob[k];
        idx2sub(k, &i, &j);
        for (kk = 0; kk < 8; kk++) {    // 借用一下minprob，界外的0当然不是雷，相当于x1
            minprob = outofbound(i, j, kk) ? 0.0 : prob[sub2idx(i + dir[kk][0], j + dir[kk][1])];
            if (minprob != I_PROB)      // I_PROB即已知格肯定不是雷，也相当于x1，和没乘一样
                temprob[k] *= (1.0 - minprob);  // 周围有雷肯定不是空白格了嘛，相当于x0
        }
    }
    for (k = 0; k < H*W; k++) {
        prob[k] -= temprob[k] * 1e-10;  // 高级最小一般为1e-5，某次中级出现1.74e-6
    }
    // 找最小概率
    for (k = mins = 0, minprob = I_PROB; k < H*W; k++) {
        if (prob[k] <= minprob) {
            if (prob[k] < minprob)  mins = 0;
            minprob = prob[k];
            min_loc[mins++] = k;
        }
    }
    if (minprob > 0.0 && minprob < 1e-6)  printf("minprob = %le\n", minprob);
    // 发现新雷先标上，有点开操作的都应放后面
    for (k = 0; k < H*W; k++) {
        if (prob[k] == 1.0 && trim[k] >= 0) {
            found_mine++;  // 不要elseif，此时可能minprob==1
            //output(trim);
            trim[k] = -9;
            idx2sub(k, &i, &j);  // unmark后trim中已没有INDEP以上的数了
            clik = mine_trim(board, trim, i, j);
            if (clik < 0)  show_error("wrong guess mine.");
        }
    }
    // 要点开了
    if (minprob == 1.0) {
        k = 0;
        mins = -1;          // 只有雷就不要点啊
        if (found_mine < MINE)  show_error("no guess, dead loop.");
    }
    else if (minprob > 0) {
        k = rand() % mins;  // 随机选一个点
        mins = k + 1;       // 只循环一次
    }
    else {
        k = 0;              // minprob<=0全点开，负数大概是减去空白格概率弄出来的
    }
    // 点！
    for (; k < mins; k++) {
        idx2sub(min_loc[k], &i, &j);
        clik = click(board, trim, i, j);
        if (clik < 0) {
            //printf("wrong guess at safe prob %6.3f%%.\n", 100 * (1 - minprob));
            //show_error("wrong guess with prob.\n");
            return -1;
        }
    }

    return found_mine;
}

// 还真不能直接融入，status用到的地方挺多
void checkSurroundTiles(int *board, int i, int j, int *status) {
    for (int k = 0; k < 8; k++)
        if (outofbound(i, j, k))
            status[k] = 0;  // 不宜用(!out)*board[]，因为可能索引已经越界，还读可能会被杀掉
        else
            status[k] = board[sub2idx(i+dir[k][0], j+dir[k][1])];
}

void markSurroundTiles(int *board, int i, int j, int *status, int markas) {
    for (int k = 0; k < 8; k++)
        if (status[k])  // outofbound的都是0
            board[sub2idx(i+dir[k][0], j+dir[k][1])] = markas;
}

void markSharedTiles(int *board, int i, int j) {
    int status[8];
    checkSurroundTiles(board, i, j, status);

    for (int k = 0; k < 8; k++)
        if (status[k] >= CLICK && status[k] < SHARE)  // SHARE区间就不用加了
            board[sub2idx(i+dir[k][0], j+dir[k][1])] += 10;
}

void unmarkSharedTiles(int *board, int i, int j) {
    int status[8];
    checkSurroundTiles(board, i, j, status);

    for (int k = 0; k < 8; k++)
        if (status[k] >= INDEP && status[k] < SHARE + 10)  // CLICK区间就不用减了
            board[sub2idx(i+dir[k][0], j+dir[k][1])] -= 10;
}

int countif(int *data, int length, int val) {
    int count = 0;
    for (int i = 0; i < length; i++)
        count += (data[i] == val);
    return count;
}

void incif(double *setdata, const int *cdata, int length, int cval, double incval) {
#pragma ivdep
    for (int i = 0; i < length; i++)
        setdata[i] += (cdata[i] == cval) * incval;
}

void incax(double *setdata, const double *xdata, int length, double alpha) {
#pragma ivdep
    for (int i = 0; i < length; i++)
        setdata[i] += alpha * xdata[i];
}

// 速算组合数的值，正确返回C(n,0)=1，要求0<=n<=12
int cnk(int n, int k) {
    int numer = 1, denom = 1;

    if (n < 0 || n < k || k < 0)  return 0;
    if (k > n - k)  k = n - k;  // min(k, n-k)

    for (int i = 0; i < k; i++)  {
        numer *= n - i;  // C *= (n - i) / (i + 1);
        denom *= i + 1;  // 但直接用上式会出分数约掉
    }
    return (numer / denom);
}

// 计算结果很大的组合数，15位有效数字，C(0,0)=1
double Cnk(int n, int k) {
    double res = 1;
    if (n < 0 || n < k || k < 0)  return 0;
    if (k > n - k)  k = n - k;  // min(k, n-k)

    for (int i = 1; i <= k; i++)
        res *= (n-i+1) / (double)i; // matlab为 (n-k+i)/i
    return res;
}

// 给出每种有意义的组合的内容，要求n<=8，k=0不应该进来
// 最终用res中的len行k列，res最后多给了一行临时空间
// void combnk(int n, int k, int (*res)[8], int len, int *c) {
//     if (len <= 0)  return;
//     if (n < 0 || n < k || k < 0)  return;  // n=0也是可以用来填的
//
//     if (k == 0) {                          // 不能区分开始就是C(n,0)还是递归来的
//         for (int i = 0; i < 8; i++)        // K之后都是乱码，当然用不到
//             res[*c][i] = res[len][i];      // c记录已有结果个数
//         (*c)++;
//     }
//
//     for (int i = n; i >= k; i--) {         // 下面是核心，最后一行是临时空间，换着数填
//         res[len][k - 1] = i - 1;           // 从右往左(k-1)，从大到小(i-1)地不重复地填0~n-1
//         combnk(i - 1, k - 1, res, len, c); // 每次把换数范围(n)缩小，k比n或i都先减到0
//     }
// }

// 要求a初始为升序地从0开始依次加一；返回1表示还有下一个
int next_combination(int n, int k, int *a) {
    int i;
    if (n <= 0 || k <= 0 || n < k)  return 0;  // n<k会死循环
    for (i = k - 1; i >= 0 && a[i] == n - k + i; i--);  // 每个位置有其最大值
    if (i < 0)  return 0;  // 排完了
    for (a[i++]++; i < k; i++)  a[i] = a[i - 1] + 1;
    return 1;
}


// 用double因为情况数可能超过uint64的1e19，11格最多已有C(8,4)^11>1e20
// 但double在2^52以上+1就没有变化了，因为IEEE用52位表示小数部分
int merge_case(double **freq_board, double *new_freq_board,
               double **case_count, double new_case_count,
               int **num_mine,      int new_num_mine,
               int *res_length, int *res_capacity, int **res_hash)
{
    if (new_case_count <= 0)  return 0;
    //if (new_case_count > 1e19)  printf("Wow, so big freq.\n");  // 一亿局真的有一次，还是XP

    int idx = (*res_hash)[new_num_mine];  // hash初始-1，new_num_mine也许越界？

    if (idx >= 0) {      // 合并到已有位置，初始化时已分配了1个
        (*case_count)[idx] += new_case_count;  // 寻址[]优先级高于取内容*
        incax(*freq_board + H*W*idx, new_freq_board, H*W, 1.0);  // freq(:,:,idx) += new_...
    }

    else if (idx < 0) {  // 放入末尾再往后一个的位置
        idx = *res_length - (**case_count == 0);  // length最初为1，第一次令idx=0以覆盖case_count[0]=0

        if (idx >= *res_capacity) {
            *res_capacity += *res_capacity / 2;  // 1.5x growth，因此capacity初值至少2
            if (idx >= *res_capacity)  *res_capacity = idx + 1;  // at least bigger than idx

            *freq_board = (double*)realloc(*freq_board, sizeof(double)**res_capacity*H*W);
            *case_count = (double*)realloc(*case_count, sizeof(double)**res_capacity);
            *num_mine   =    (int*)realloc(*num_mine,   sizeof(int)**res_capacity);
            if (!*freq_board || !*case_count || !*num_mine)  return -1;
        }

        *res_length = idx + 1;  // 放在realloc后不然读取越界啊
        (*res_hash)[new_num_mine] = idx;

        (*num_mine)[idx] = new_num_mine;
        (*case_count)[idx] = new_case_count;  // freq_board(:,:,idx) = new_freq_board;
        memcpy(*freq_board + H*W*idx, new_freq_board, sizeof(double)*H*W);
    }

    return 0;
}

// 非公区对同一个雷数各格概率都一样(记下每格概率)，但由于相互独立，[公共区的一种排布]可搭配[非公区各自排出所有可能]
// 公共区需要试验该雷数下每种排布(记下区内每格是雷不是)并接着验证下一个数字
// 直到验证完所有数字时，非公区用概率乘上[非公区各自排出所有可能的总数]，公共区的次数就是每格有没有雷乘总数
// 这样就记下了一个雷数下的其中一类("总数"种情况)中的每格是雷的[次数]，进行merge_case
// 另外，只有depth、cur_mine、freq_coef会每层不同；inc_freq_board、freq_board、num_mine仅最后一层用到
int try_mine_recur(int *board, int *ii, int *jj, int max_mine, int max_depth,
                   int depth, int cur_mine, double freq_coef, double *freq_unshared,
                   double *inc_freq_board, res_info *info)
{
    if (depth >= max_depth) {
        memset(inc_freq_board, 0, sizeof(double)*H*W);
        incif(inc_freq_board, board, H*W, H_MINE, freq_coef);
        incax(inc_freq_board, freq_unshared, H*W, freq_coef);

        // ->优先级高于&
        return merge_case(&info->freq_board, inc_freq_board,
                          &info->case_count, freq_coef,
                          &info->num_mine,   cur_mine,
                          &info->res_length, &info->res_capacity, &info->res_hash);
    }

    int i = ii[depth], j = jj[depth], need_mine, retVal, status[8];
    int num_shared = 0, shared_idx[8], shared_dir[8], comb_shared[8];
    int num_unshared = 0, unshared_mine, cnk_unshared;
    double unshared_prob;


    checkSurroundTiles(board, i, j, status);
    need_mine = -board[sub2idx(i, j)] - countif(status, 8, H_MINE);  // board已trim过，已知为负数

    // 便于向量化的放一起
    for (int k = 0; k < 8; k++) {
        shared_dir[k] = (ISSHARE(status[k]));  // true=1
        num_unshared += (ISINDEP(status[k]));
    }

    for (int k = 0; k < 8; k++) {
        if (shared_dir[k])
            shared_idx[num_shared++] = k;
    }

    // 命中率高的条件放前面；要放的雷数超过格子数或已有雷数超过总雷数
    if (need_mine > num_shared + num_unshared || need_mine < 0)  return 0;
    if (cur_mine + need_mine > max_mine || cur_mine > max_mine)  return 0;

    // 留shared_mine个雷去放公共区，剩下雷的先在非公区独立地排排算算概率
    for (int shared_mine = 0; shared_mine <= need_mine; shared_mine++) {

        cnk_unshared = 0;
        unshared_mine = 0;

        // 计算非公共区格概率，跳过使C(n,k)无意义的数值(格子要够)
        if (num_unshared >= need_mine - shared_mine) {
            unshared_mine = need_mine - shared_mine;
            cnk_unshared = cnk(num_unshared, unshared_mine);

            // 有格有雷，cnk>0算频数和概率；有格没雷，非公区就全不放雷
            unshared_prob = cnk(num_unshared-1, unshared_mine-1) / (double)cnk_unshared;

            for (int k = 0; k < 8; k++)  // 不在非公区的格就0，和背景一样
                if (ISINDEP(status[k]))  // unshared
                    freq_unshared[sub2idx(i+dir[k][0], j+dir[k][1])] = unshared_prob;
        }

        if (cnk_unshared == 0)  continue;  // 0说明非公区格子不够，连上面的if都没进，不满足此数字

        // 计算公共区格概率，跳过无意义数值
        if (num_shared >= shared_mine) {
            for (int k = 0; k < 8; k++)  comb_shared[k] = k;  // 准备从头排列雷

            // 把周围格可能的放雷方式列出来，一个个试
            do {
                markSurroundTiles(board, i, j, shared_dir, H_SAFE);

                // shared_mine为0不放雷，但也得继续递归
                for (int m, k = 0; k < shared_mine; k++) {
                    m = shared_idx[comb_shared[k]];  // 放的雷会覆盖shared_dir的一些点
                    board[sub2idx(i+dir[m][0], j+dir[m][1])] = H_MINE;
                }

                retVal =
                try_mine_recur(board, ii, jj, max_mine, max_depth,
                               depth+1, cur_mine+unshared_mine+shared_mine, freq_coef*cnk_unshared,
                               freq_unshared, inc_freq_board, info);
                if (retVal != 0)  return -1;

                // 恢复现场，不能像以前之间设为SHARED，现在trim还留着周围已知雷数信息
                for (int k = 0; k < 8; k++)
                    if (shared_dir[k])
                        board[sub2idx(i+dir[k][0], j+dir[k][1])] = status[k];

            } while (next_combination(num_shared, shared_mine, comb_shared));
        }

        // 恢复非公共区格的现场...注意忘记判越界会heap corruption
        if (num_unshared >= need_mine - shared_mine) {  // 不等号方向啊前面if改了这里忘了啊
            for (int k = 0; k < 8; k++)  // 非公区必互不重合，退出时直接周围全0
                if (ISINDEP(status[k]))  // 不然会越界，这里没有两圈墙啊
                    freq_unshared[sub2idx(i+dir[k][0], j+dir[k][1])] = 0.0;
        }
    }

    return 0;
}
