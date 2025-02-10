// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <omp.h> // 引入 OpenMP 头文件

using namespace Rcpp;
using namespace std;

// 定义分布类型枚举
enum DistributionType { BINOMIAL, POISSON, NEG_BINOMIAL };

// 结构体定义，用于存储候选模型的参数和评估标准
struct CandidateModel {
    double criterion; // 用于存储AIC/BIC/AICc值
    long pos_start;    // 序列的起始位置
    long pos_end;      // 序列的结束位置
    long cs;           // 聚类的起始位置
    long ce;           // 聚类的结束位置
    double p0;         // 非聚类区域的死亡率
    double pc;         // 聚类区域的死亡率
    double InL0;       // Null model 的对数似然
    double InL;        // 当前模型的对数似然
    double AIC0;       // Null model 的 AIC
    double AIC;        // 当前模型的 AIC
    double AICc0;      // Null model 的 AICc
    double AICc;       // 当前模型的 AICc
    double BIC0;       // Null model 的 BIC
    double BIC;        // 当前模型的 BIC
    double weight;     // 模型平均化的权重

    CandidateModel() {}
    CandidateModel(long ps, long pe, long c_s, long c_e, double p0_, double pc_,
                   double InL0_, double InL_, double AIC0_, double AIC_,
                   double AICc0_, double AICc_, double BIC0_, double BIC_) :
        pos_start(ps), pos_end(pe), cs(c_s), ce(c_e), p0(p0_), pc(pc_),
        InL0(InL0_), InL(InL_), AIC0(AIC0_), AIC(AIC_), 
        AICc0(AICc0_), AICc(AICc_), BIC0(BIC0_), BIC(BIC_), weight(0.0) {}
};

// 计算 Binomial 分布的对数似然
double binomial_loglik(long deaths, long pop, double p) {
    if (p <= 0.0 || p >= 1.0) {
        return -std::numeric_limits<double>::infinity();
    }
    return lgamma(pop + 1) - lgamma(deaths + 1) - lgamma(pop - deaths + 1) +
           deaths * log(p) + (pop - deaths) * log(1.0 - p);
}

// 计算 Null Model 的对数似然和信息准则
void calc_null_model(const std::vector<long> &deaths, const std::vector<long> &pop,
                    long pos_start, long pos_end,
                    double &InL0, double &AIC0, double &AICc0, double &BIC0,
                    double &p_all) {
    long total_deaths = 0;
    long total_pop = 0;
    for(long i = pos_start; i <= pos_end; i++) {
        total_deaths += deaths[i];
        total_pop += pop[i];
    }
    if(total_pop == 0){
        p_all = 0.0;
    }
    else{
        p_all = static_cast<double>(total_deaths) / static_cast<double>(total_pop);
    }

    // 计算 Null Model 的对数似然
    InL0 = 0.0;
    for(long i = pos_start; i <= pos_end; i++) {
        InL0 += binomial_loglik(deaths[i], pop[i], p_all);
    }

    // 计算 AIC, AICc, BIC
    int k = 1; // p_all 是一个参数
    double l = static_cast<double>(pos_end - pos_start + 1);
    AIC0 = -2.0 * InL0 + 2.0 * k;
    if(l - k - 1 > 0){
        AICc0 = AIC0 + 2.0 * k * (k + 1) / (l - k - 1);
    }
    else{
        AICc0 = AIC0;
    }
    BIC0 = -2.0 * InL0 + k * log(l);
}

// 计算 Poisson 分布的对数似然
double poisson_loglik(long deaths, double lambda) {
    if(lambda <= 0.0){
        return -std::numeric_limits<double>::infinity();
    }
    return deaths * log(lambda) - lambda - lgamma(deaths + 1);
}

// 计算 Negative Binomial 分布的对数似然
double neg_binomial_loglik(long deaths, double p, double r) {
    if(p <= 0.0 || p >= 1.0 || r <= 0.0){
        return -std::numeric_limits<double>::infinity();
    }
    return lgamma(deaths + r) - lgamma(r) - lgamma(deaths + 1) +
           r * log(1.0 - p) + deaths * log(p);
}

// 通用的对数似然计算函数，根据分布类型选择相应的计算方法
double calc_loglik(const std::vector<long> &deaths, const std::vector<long> &pop, 
                  long pos_start, long pos_end, long cs, long ce, 
                  double p0, double pc, DistributionType dist_type) {
    double ll = 0.0;
    for(long i = pos_start; i <= pos_end; i++) {
        double curp = (i >= cs && i <= ce) ? pc : p0;
        if(dist_type == BINOMIAL){
            ll += binomial_loglik(deaths[i], pop[i], curp);
        }
        else if(dist_type == POISSON){
            ll += poisson_loglik(deaths[i], curp);
        }
        else if(dist_type == NEG_BINOMIAL){
            double r = 1.0; // 过度参数，可以调整或优化
            ll += neg_binomial_loglik(deaths[i], curp, r);
        }
    }
    return ll;
}

// 寻找给定子序列的最佳模型
CandidateModel find_best_model(const std::vector<long> &deaths, const std::vector<long> &pop, 
                               long pos_start, long pos_end, int criterion_type, 
                               DistributionType dist_type){
    // 计算 Null Model
    double InL0, AIC0, AICc0, BIC0, p_all;
    calc_null_model(deaths, pop, pos_start, pos_end, InL0, AIC0, AICc0, BIC0, p_all);

    // 初始化最佳模型为 Null Model
    double best_criterion = std::numeric_limits<double>::infinity();
    CandidateModel best_model(pos_start, pos_end, pos_start, pos_end, p_all, p_all, 
                              InL0, InL0, AIC0, AIC0, AICc0, AICc0, BIC0, BIC0);

    // 并行遍历所有 (cs, ce) 组合
    #pragma omp parallel
    {
        double local_best_criterion = std::numeric_limits<double>::infinity();
        CandidateModel local_best_model;

        #pragma omp for nowait
        for(long cs = pos_start; cs <= pos_end; cs++) {
            for(long ce = cs + 1; ce <= pos_end; ce++) {
                // 计算分段参数
                long deaths_s = 0, pop_s = 0;
                for(long i = pos_start; i < cs; i++) { deaths_s += deaths[i]; pop_s += pop[i]; }
                long deaths_c = 0, pop_c = 0;
                for(long i = cs; i <= ce; i++) { deaths_c += deaths[i]; pop_c += pop[i]; }
                long deaths_e = 0, pop_e = 0;
                for(long i = ce + 1; i <= pos_end; i++) { deaths_e += deaths[i]; pop_e += pop[i]; }

                double p0_est, pc_est;
                long pop_out = pop_s + pop_e;
                long deaths_out = deaths_s + deaths_e;
                if(pop_out == 0){
                    p0_est = p_all;
                }
                else{
                    p0_est = static_cast<double>(deaths_out) / static_cast<double>(pop_out);
                }

                // 根据分布类型计算 p_in
                if(dist_type == BINOMIAL){
                    if(pop_c == 0){
                        pc_est = p_all;
                    }
                    else{
                        pc_est = static_cast<double>(deaths_c) / static_cast<double>(pop_c);
                    }
                }
                else {
                    // 对于 Poisson 和 Negative Binomial，p0 和 pc 代表不同参数
                    pc_est = static_cast<double>(deaths_c) / static_cast<double>(pop_c);
                }

                // 计算当前模型的对数似然
                double InL = calc_loglik(deaths, pop, pos_start, pos_end, cs, ce, p0_est, pc_est, dist_type);

                // 计算信息准则
                double k = 2.0; // p0 和 pc 是两个参数
                double l = static_cast<double>(pos_end - pos_start + 1);
                double AIC = -2.0 * InL + 2.0 * k;
                double AICc = AIC;
                if(l - k - 1 > 0){
                    AICc += 2.0 * k * (k + 1) / (l - k - 1);
                }
                double BIC = -2.0 * InL + k * log(l);

                double cri;
                if(criterion_type == 0){ // BIC
                    cri = BIC;
                }
                else if(criterion_type == 1){ // AIC
                    cri = AIC;
                }
                else { // AICc
                    cri = AICc;
                }

                if(cri < local_best_criterion){
                    local_best_criterion = cri;
                    local_best_model = CandidateModel(pos_start, pos_end, cs, ce, p0_est, pc_est, 
                                                      InL0, InL, AIC0, AIC, AICc0, AICc, BIC0, BIC);
                }
            }
        }

        // 保护性地更新全局最佳模型
        #pragma omp critical
        {
            if(local_best_criterion < best_criterion){
                best_criterion = local_best_criterion;
                best_model = local_best_model;
            }
        }
    }

    return best_model;
}

// 递归的 divide-and-conquer 过程，根据最佳分割点继续细分
void cluster_subseq(const std::vector<long> &deaths, const std::vector<long> &pop, 
                   long pos_start, long pos_end, int criterion_type, 
                   DistributionType dist_type,
                   std::vector<CandidateModel> &selected_models){

    // 寻找本段序列最佳模型
    CandidateModel best_model = find_best_model(deaths, pop, pos_start, pos_end, criterion_type, dist_type);

    // 判断是否有显著聚类
    // 简化判断条件：如果 best_model 和 null model 相同，则无聚类
    // 可以根据实际需要添加更复杂的检验
    if(best_model.cs == pos_start && best_model.ce == pos_end && fabs(best_model.p0 - best_model.pc) < 1e-10){
        // 无显著分段
        selected_models.push_back(best_model);
        return;
    }

    // 有显著分段则记录该 model，并递归左右和中间段
    selected_models.push_back(best_model);

    // 左侧区间
    if(best_model.cs > pos_start + 1){
        cluster_subseq(deaths, pop, pos_start, best_model.cs - 1, criterion_type, dist_type, selected_models);
    }

    // 中间区间
    if(best_model.ce > best_model.cs){
        cluster_subseq(deaths, pop, best_model.cs, best_model.ce, criterion_type, dist_type, selected_models);
    }

    // 右侧区间
    if(best_model.ce < pos_end - 1){
        cluster_subseq(deaths, pop, best_model.ce + 1, pos_end, criterion_type, dist_type, selected_models);
    }

}

// [[Rcpp::export]]
List runClusterSurvival(IntegerVector deaths_r, IntegerVector pop_r, int criterion_type=0, int dist_type_code=0){
    // criterion_type: 0=BIC,1=AIC,2=AICc
    // dist_type_code: 0=BINOMIAL,1=POISSON,2=NEG_BINOMIAL

    // 转换输入数据
    std::vector<long> deaths(deaths_r.begin(), deaths_r.end());
    std::vector<long> pop(pop_r.begin(), pop_r.end());

    // 确定分布类型
    DistributionType dist_type;
    if(dist_type_code == 0){
        dist_type = BINOMIAL;
    }
    else if(dist_type_code == 1){
        dist_type = POISSON;
    }
    else{
        dist_type = NEG_BINOMIAL;
    }

    // 存储选定的模型
    std::vector<CandidateModel> selected_models;

    // 执行分段聚类
    cluster_subseq(deaths, pop, 0, static_cast<long>(deaths.size()) - 1, criterion_type, dist_type, selected_models);

    // 计算模型平均化
    // 计算最小信息准则值
    double min_criterion = std::numeric_limits<double>::infinity();
    for(auto &model : selected_models){
        if(criterion_type == 0){ // BIC
            if(model.BIC < min_criterion){
                min_criterion = model.BIC;
            }
        }
        else if(criterion_type == 1){ // AIC
            if(model.AIC < min_criterion){
                min_criterion = model.AIC;
            }
        }
        else { // AICc
            if(model.AICc < min_criterion){
                min_criterion = model.AICc;
            }
        }
    }

    // 计算模型权重
    double sum_weights = 0.0;
    for(auto &model : selected_models){
        double delta;
        if(criterion_type == 0){
            delta = model.BIC - min_criterion;
        }
        else if(criterion_type == 1){
            delta = model.AIC - min_criterion;
        }
        else{
            delta = model.AICc - min_criterion;
        }
        model.weight = exp(-0.5 * delta);
        sum_weights += model.weight;
    }

    // 归一化权重
    for(auto &model : selected_models){
        model.weight /= sum_weights;
    }

    // 计算加权平均 p0 和 pc
    double avg_p0 = 0.0, avg_pc = 0.0;
    for(auto &model : selected_models){
        avg_p0 += model.p0 * model.weight;
        avg_pc += model.pc * model.weight;
    }

    // 生成每个位置的 p(i)
    int N = deaths.size();
    std::vector<double> p_i(N, 0.0);
    for(auto &model : selected_models){
        for(int i = model.pos_start; i <= model.pos_end; i++){
            if(i >= model.cs && i <= model.ce){
                p_i[i] += model.pc * model.weight;
            }
            else{
                p_i[i] += model.p0 * model.weight;
            }
        }
    }

    // 将 p_i 转换为 R 的 NumericVector 并添加到结果中
    NumericVector p_i_r(N);
    for(int i = 0; i < N; i++){
        p_i_r[i] = p_i[i];
    }

    // 将分段点及其对应参数提取到 R 的向量中
    IntegerVector cs_vec, ce_vec;
    NumericVector p0_vec, pc_vec;
    for(auto &model : selected_models){
        cs_vec.push_back(model.cs + 1); // 转为1-based
        ce_vec.push_back(model.ce + 1);
        p0_vec.push_back(model.p0);
        pc_vec.push_back(model.pc);
    }

    // 打包结果
    List result;
    result["cs"] = cs_vec;
    result["ce"] = ce_vec;
    result["p0"] = p0_vec;
    result["pc"] = pc_vec;
    result["avg_p0"] = avg_p0;
    result["avg_pc"] = avg_pc;
    result["p_i"] = p_i_r;

    return result;
}
