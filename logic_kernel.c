#include "ms.h"

extern int dir[][2];
static int dir2[][2] = { { -2, -2 }, { -2, -1 }, { -2, 0 }, { -2, 1 }, { -2, 2 },   //  1  2  3  4  5
                         { -1, -2 }, { -1, -1 }, { -1, 0 }, { -1, 1 }, { -1, 2 },   //  6  7  8  9 10
                         {  0, -2 }, {  0, -1 },            {  0, 1 }, {  0, 2 },   // 11 12    13 14
                         {  1, -2 }, {  1, -1 }, {  1, 0 }, {  1, 1 }, {  1, 2 },   // 15 16 17 18 19
                         {  2, -2 }, {  2, -1 }, {  2, 0 }, {  2, 1 }, {  2, 2 } }; // 20 21 22 23 24

int outofbound2(int i, int j, int k)
{ return (i + dir2[k][0] < 0 || i + dir2[k][0] >= H || j + dir2[k][1] < 0 || j + dir2[k][1] >= W); }


int logic_kernel(int *board, int *trim, int found_mine)
{
    int clik;

    int mark_1(int*, int*, int*);
    int mark_2(int*, int*, int*);

    do {
        // 1-推理，mark1推理在一次推出得多的时候可连用几次
        do {
            clik = mark_1(board, trim, &found_mine);
            if (clik < 0)  return -1;
        } while (clik > 0);
        //output(trim);

        // 2-推理，这个比较耗时所以少跑一点，多推点简单的发现新线索
        clik = mark_2(board, trim, &found_mine);  // 现在发现雷时mark1/2机制是一样的
        if (clik < 0)  return -1;                 // 不用像matlab最后再来一次mark1
    } while (clik > 0);

    return found_mine;
}


// 把指定位置的雷周围数字值减1【未知格先标记下这个以后点开时数字要减1】减到0顺便点开那个数字的周围
// PS 点开格子(也会改数字值)不能打断此处改值过程，不然会漏减或多减新点开的数字，遂把可点开的先入队
//    因此mark2前和guess前要出现一次[没开新格的trim]，即调用此函数的mark1/2里对此处返回值累加得0
int mine_trim(int *board, int *trim, int row, int col)
{
    int k, r, s, t, ii, jj, kk, clik, area, len, pos[64];

    len = 0;
    area = 0;

    for (k = 0; k < 8; k++) {  // 雷周围一圈数字值减一
        if (outofbound(row, col, k))  continue;
        r = sub2idx(row + dir[k][0], col + dir[k][1]);

        t = -trim[r];
        if (!(t > 0 && t < 9 || -t >= CLICK))  continue;

        t = -(++trim[r]);      // 已知数字是负的；可点格则恰需增加
        if (t == 0) {          // 减成空白格要特殊处理
            trim[r] = -10;     // 重设数字+点开周围

            ii = row + dir[k][0];
            jj = col + dir[k][1];
            for (kk = 0; kk < 8; kk++) {  // 点开周围，注意两重位移
                if (outofbound(ii, jj, kk))  continue;
                s = sub2idx(ii + dir[kk][0], jj + dir[kk][1]);

                if (trim[s] >= 0) {       // >=0表示未知，重复入队可达16格(雷旁一圈1，再外圈全可点时)
                    pos[len++] = sub2idx(ii + dir[kk][0], jj + dir[kk][1]);  // 虽然重复也不会出错
                }
            }
        }
    }

    for (k = 0; k < len; k++) {
        idx2sub(pos[k], &ii, &jj);
        clik = click(board, trim, ii, jj);
        if (clik < 0)  return -1;
        area += clik;
    }

    return area;
}

// 根据1格数字推断旁边的雷和安全格，碰到雷顺便trim
int mark_1(int *board, int *trim, int *mines)
{
    int i, j, k, r, t, clik, area;

    int count_unknown(int*, int, int);

    area = 0;  // 1-推断可重新找出所有雷，但要看外面历史进程的要求

    for (i = 0; i < H; i++) {
        for (j = 0; j < W; j++) {
            t = -trim[sub2idx(i, j)];
            if (!(t > 0 && t < 9))  continue;  // 已知数字才判断

            if (t == count_unknown(trim, i, j)) {  // 该标雷
                for (k = 0; k < 8; k++) {
                    if (outofbound(i, j, k))  continue;
                    r = sub2idx(i + dir[k][0], j + dir[k][1]);

                    if (trim[r] >= 0) {  // 即未知格子(10)
                        trim[r] = -9;    // 雷被发现后应mine_trim，使雷旁边的数字变小
                        (*mines)++;      // 未知格数与数字同降，count_unknown不会错
                        clik = mine_trim(board, trim, i + dir[k][0], j + dir[k][1]);
                        if (clik < 0)  show_error("wrong mark_1.");
                        area += clik;
                    }
                }
            }
        }
    }

    return area;
}

// 根据2格数字推断旁边的雷和安全格，碰到雷顺便trim
int mark_2(int *board, int *trim, int *mines)
{
    int k, r, s, clik, area;
    int i1, i2, j1, j2, num1, num2, rest1, rest2, status1[8], status2[8];
    int shared_tile, shared_mine_max, shared_mine_min;

    area = 0;  // mines是在mark1基础上再找，不必归零

    for (i1 = 0; i1 < H; i1++) {
        for (j1 = 0; j1 < W; j1++) {
            for (k = 0; k < 24; k++) {  // 找中心数字的两圈邻域

                num1 = -trim[sub2idx(i1, j1)];  // 标雷后包括中心数字都会变！不更新标雷错啊
                if (!(num1 > 0 && num1 < 9))  break;  // 中心数字没了就跳出k循环换个数字

                if (outofbound2(i1, j1, k))  continue;
                i2 = i1 + dir2[k][0];
                j2 = j1 + dir2[k][1];
                num2 = -trim[sub2idx(i2, j2)];
                if (!(num2 > 0 && num2 < 9))  continue;  // 俩都得是数字

                markSharedTiles(trim, i1, j1);   // 在循环外unmark时会多减
                markSharedTiles(trim, i2, j2);   // i1还是i2? 复制来要改啊
                checkSurroundTiles(trim, i1, j1, status1);
                checkSurroundTiles(trim, i2, j2, status2);
                unmarkSharedTiles(trim, i2, j2); // 搞完就取消，不然出错
                unmarkSharedTiles(trim, i1, j1); // 有INDEP这些容易错

                rest1 = rest2 = shared_tile = 0;
                for (s = 0; s < 8; s++) {
                    rest1 += ISINDEP(status1[s]);
                    rest2 += ISINDEP(status2[s]);
                    shared_tile += ISSHARE(status1[s]);  // 俩shared数量显然一样
                }
                if (rest1 == 0 || shared_tile == 0)  continue;  // 下面for只弄rest1的格子

                shared_mine_max = min(min(num1, num2), shared_tile);
                shared_mine_min = max(max(num1 - rest1, num2 - rest2), 0);

                if (shared_mine_min >= num1) {  // 数字的剩余区无雷时，公共区需要雷数最多
                    for (s = 0; s < 8; s++) {   // 若公共区最少雷数已满足此需求，则剩余区都安全
                        if (ISINDEP(status1[s])) {  // unmark后trim中已没有INDEP以上的数了
                            clik = click(board, trim, i1 + dir[s][0], j1 + dir[s][1]);
                            if (clik < 0)  show_error("wrong mark_2 safe.");
                            area += clik;
                        }
                    }
                }

                if (shared_mine_max <= num1 - rest1) {  // 实际不会有小于或max取num1(会恰好=num2)
                    for (s = 0; s < 8; s++) {       // <=右边假设剩余区都是雷，此时公共区要放的雷最少
                        if (ISINDEP(status1[s])) {  // 若公共区最大雷数还不满足需要，则剩余区只能全雷
                            r = sub2idx(i1 + dir[s][0], j1 + dir[s][1]);  // k还是s复制来要改啊

                            if (trim[r] >= 0) {  // 防止重复标雷
                                trim[r] = -9;    // mark2找到雷也要trim+点开
                                (*mines)++;      // 不trim的话给guess的数字就不对了啊
                                clik = mine_trim(board, trim, i1 + dir[s][0], j1 + dir[s][1]);
                                if (clik < 0)  show_error("wrong mark_2 mine.");
                                area += clik;
                            }
                        }
                    }
                }
            }
        }
    }

    return area;
}


int count_unknown(int *board, int row, int col)
{
    int k, r, unknown;
    for (k = unknown = 0; k < 8; k++) {
        if (outofbound(row, col, k))  continue;
        r = sub2idx(row + dir[k][0], col + dir[k][1]);
        unknown += (board[r] >= 0);  // >=0表示未知格
    }
    return unknown;
}
