#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// 设置最小分段长度（观察点数）和注入的惩罚参数 lambda（可由用户输入） 
const long MIN_SEGMENT_LENGTH = 3;

// 二项分布对数似然函数
// deaths ~ Binomial(total, p)
// 对数似然： l = deaths*log(p) + (total-deaths)*log(1-p)
double binomial_loglik(long deaths, long total, double p) {
    if(p <= 0.0 || p >= 1.0) return -std::numeric_limits<double>::infinity();
    return deaths * log(p) + (total - deaths) * log(1.0 - p);
}

// 计算区间 [i, j] 内成本，成本定义为 -2*loglikelihood
// 对于该区间，采用 MLE p = (sum_{k=i}^j deaths[k])/(sum_{k=i}^j pop[k])
double segment_cost(int i, int j, const vector<long>& deaths, const vector<long>& pop) {
    long totd = 0, totpop = 0;
    for (int k = i; k <= j; k++) {
        totd += deaths[k];
        totpop += pop[k];
    }
    if(totpop == 0) return 0.0;
    double p = (double) totd / totpop;
    double ll = 0.0;
    for (int k = i; k <= j; k++) {
        ll += binomial_loglik(deaths[k], pop[k], p);
    }
    return -2.0 * ll;
}

// 动态规划分段函数，采用全局优化，目标为最小化总体成本 + λ * (#段)
// 返回分段变化点（以 0 开始的索引）以及每个分段的 MLE p 值。
//
// 参数说明：
//   deaths_r: 每个时间点的死亡事件数
//   pop_r: 每个时间点的总体风险数
//   lambda: 每个分段引入的惩罚（可以理解为平滑参数，较大 lambda 越不容易分段）
//
// [[Rcpp::export]]
List runClusterSurvivalDP(IntegerVector deaths_r, IntegerVector pop_r, double lambda) {
    int n = deaths_r.size();
    vector<long> deaths(n), pop(n);
    for (int i = 0; i < n; i++) {
        deaths[i] = deaths_r[i];
        pop[i] = pop_r[i];
    }
    
    // 动态规划数组 D[i]: 表示前 i 个观测的最小累计成本
    vector<double> D(n + 1, std::numeric_limits<double>::infinity());
    // last_cut[i] 记录 i 处最优分段的上一个切分点
    vector<int> last_cut(n + 1, -1);
    
    D[0] = 0;
    
    // DP求解：对于每个终点 j，从 1 到 n
    for (int j = 1; j <= n; j++) {
        // i 表示上一个切分点，区间 [i, j-1]为当前段
        for (int i = 0; i <= j - MIN_SEGMENT_LENGTH; i++) {
            double cost = segment_cost(i, j - 1, deaths, pop);
            double total = D[i] + cost + lambda;
            if (total < D[j]) {
                D[j] = total;
                last_cut[j] = i;
            }
        }
    }
    
    // 反向追踪获得分段变化点
    vector<int> change_points;
    int pos = n;
    while (pos > 0) {
        int cp = last_cut[pos];
        if (cp < 0) break;
        change_points.push_back(cp);
        pos = cp;
    }
    reverse(change_points.begin(), change_points.end());
    
    // 根据变化点确定每个分段的起始和结束位置以及 MLE p 值
    vector<int> seg_start;
    vector<int> seg_end;
    vector<double> seg_p;
    
    int start = 0;
    for (size_t i = 0; i < change_points.size(); i++) {
        int end = change_points[i] - 1;
        if (end < start) continue;
        seg_start.push_back(start);
        seg_end.push_back(end);
        long totd = 0, totpop = 0;
        for (int k = start; k <= end; k++){
            totd += deaths[k];
            totpop += pop[k];
        }
        double p = (totpop == 0) ? 0.0 : (double) totd / totpop;
        seg_p.push_back(p);
        start = change_points[i];
    }
    if (start < n) {
        seg_start.push_back(start);
        seg_end.push_back(n - 1);
        long totd = 0, totpop = 0;
        for (int k = start; k < n; k++){
            totd += deaths[k];
            totpop += pop[k];
        }
        double p = (totpop == 0) ? 0.0 : (double) totd / totpop;
        seg_p.push_back(p);
    }
    
    // 将结果打包返回
    IntegerVector change_points_r(change_points.begin(), change_points.end());
    IntegerVector seg_start_r(seg_start.begin(), seg_start.end());
    IntegerVector seg_end_r(seg_end.begin(), seg_end.end());
    NumericVector seg_p_r(seg_p.begin(), seg_p.end());
    
    return List::create(Named("change_points") = change_points_r,
                        Named("segment_start") = seg_start_r,
                        Named("segment_end") = seg_end_r,
                        Named("segment_p") = seg_p_r,
                        Named("total_cost") = D[n]);
}
