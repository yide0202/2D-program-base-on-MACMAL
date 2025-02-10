// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <omp.h>

using namespace Rcpp;
using namespace std;

/*
--------------------------------------------------------------------------------
(一) Binomial+Logistic 梯度与优化核心
--------------------------------------------------------------------------------
*/

// 1. logistic 函数
inline double logistic(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

// 2. Binomial + logistic 对数似然
inline double binomial_loglik(long deaths, long pop, double p){
    if(p <= 0.0 || p >= 1.0){
        return -std::numeric_limits<double>::infinity();
    }
    // 忽略组合项常量，对梯度无影响
    return deaths * std::log(p) + (pop - deaths)*std::log(1.0 - p);
}

// 3. calc_logit: 给定 beta(K+1) + xrow(K), 计算 z = beta0 + Σ( beta_k * x_k )
inline double calc_logit(const std::vector<double> &beta, const std::vector<double> &xrow){
    double val = beta[0]; // 截距
    for(size_t k=0; k < xrow.size(); k++){
        val += beta[k+1]*xrow[k];
    }
    return val;
}

// 4. calcLogLikAndGrad: 对给定样本索引, 计算对数似然+梯度
void calcLogLikAndGrad(const std::vector<long> &deaths,
                       const std::vector<long> &pop,
                       const NumericMatrix &X,
                       const std::vector<long> &indices,
                       const std::vector<double> &beta,
                       double &logLik,
                       std::vector<double> &grad)
{
    logLik = 0.0;
    std::fill(grad.begin(), grad.end(), 0.0);

    for(auto i : indices){
        // xrow
        std::vector<double> xrow(X.ncol());
        for(int k=0; k < X.ncol(); k++){
            xrow[k] = X(i, k);
        }

        double z_i = calc_logit(beta, xrow);
        double p_i = logistic(z_i);

        double ll_i = binomial_loglik(deaths[i], pop[i], p_i);
        logLik += ll_i;

        double resid = (double)deaths[i] - (double)pop[i]*p_i;
        grad[0] += resid; 
        for(size_t d=0; d < xrow.size(); d++){
            grad[d+1] += resid * xrow[d];
        }
    }
}

// 5. gradientDescentBinomial: 用固定学习率梯度下降来拟合 Beta
void gradientDescentBinomial(const std::vector<long> &deaths,
                             const std::vector<long> &pop,
                             const NumericMatrix &X,
                             const std::vector<long> &indices,
                             std::vector<double> &beta,
                             int maxIter,
                             double lr)
{
    int Kplus1 = beta.size();
    std::vector<double> grad(Kplus1, 0.0);

    for(int iter=0; iter < maxIter; iter++){
        double logLik = 0.0;
        calcLogLikAndGrad(deaths, pop, X, indices, beta, logLik, grad);

        // 更新
        for(int j=0; j < Kplus1; j++){
            beta[j] += lr * grad[j];
        }
    }
}

/*
--------------------------------------------------------------------------------
(二) 分段聚类 + 多维协变量 主逻辑
--------------------------------------------------------------------------------
*/

// ModelParams: 存储分段 (cs, ce) 及回归系数等
struct ModelParams {
    long pos_start;
    long pos_end;
    long cs;
    long ce;

    // 2组(非聚类/聚类), each has (K+1)个参数
    vector<vector<double>> beta;

    double logLik;
    double AIC;
    double BIC;

    ModelParams(int K=0){
        beta.resize(2, vector<double>(K+1, 0.0));
        logLik = -std::numeric_limits<double>::infinity();
        AIC = BIC = 0.0;
    }
};

// 计算对数似然 (考虑分段区间)
double computeLogLik(const vector<long> &deaths,
                     const vector<long> &pop,
                     const NumericMatrix &X,
                     long pos_start,
                     long pos_end,
                     long cs,
                     long ce,
                     const vector<double> &beta0,
                     const vector<double> &beta1)
{
    double ll = 0.0;
    for(long i=pos_start; i <= pos_end; i++){
        vector<double> xrow(X.ncol());
        for(int k=0; k < X.ncol(); k++){
            xrow[k] = X(i, k);
        }
        double z_i = ((i >= cs && i <= ce)
                      ? calc_logit(beta1, xrow)
                      : calc_logit(beta0, xrow));
        double p_i = logistic(z_i);
        ll += binomial_loglik(deaths[i], pop[i], p_i);
    }
    return ll;
}

// estimateCoeffs: 对 "非聚类区间" 与 "聚类区间" 分别做梯度下降
void estimateCoeffs(const vector<long> &deaths,
                    const vector<long> &pop,
                    const NumericMatrix &X,
                    long pos_start, long pos_end,
                    long cs, long ce,
                    vector<double> &beta0,
                    vector<double> &beta1)
{
    // idx_nonClust, idx_clust
    vector<long> idx_nonClust, idx_clust;
    idx_nonClust.reserve(pos_end-pos_start+1);
    idx_clust.reserve(pos_end-pos_start+1);

    for(long i=pos_start; i<=pos_end; i++){
        if(i>=cs && i<=ce){
            idx_clust.push_back(i);
        } else {
            idx_nonClust.push_back(i);
        }
    }

    // 分别调用 gradientDescentBinomial
    gradientDescentBinomial(deaths, pop, X, idx_nonClust, beta0, 100, 1e-5);
    gradientDescentBinomial(deaths, pop, X, idx_clust,    beta1, 100, 1e-5);
}

// 遍历 (cs, ce) 搜索最优分段 + 回归系数
ModelParams find_best_model_with_covariates(const vector<long> &deaths,
                                            const vector<long> &pop,
                                            const NumericMatrix &X,
                                            long pos_start,
                                            long pos_end,
                                            int K,
                                            int criterion_type)
{
    double bestCrit = std::numeric_limits<double>::infinity();
    ModelParams bestM(K);
    int N = pos_end - pos_start + 1;

    #pragma omp parallel
    {
        double local_best = std::numeric_limits<double>::infinity();
        ModelParams localModel(K);

        #pragma omp for nowait
        for(long cs=pos_start; cs<=pos_end; cs++){
            for(long ce=cs+1; ce<=pos_end; ce++){
                // 初始化
                vector<double> b0(K+1, 0.0);
                vector<double> b1(K+1, 0.0);

                estimateCoeffs(deaths, pop, X, pos_start, pos_end, cs, ce,
                               b0, b1);

                double ll = computeLogLik(deaths, pop, X, pos_start, pos_end,
                                          cs, ce, b0, b1);

                int num_params = 2*(K+1);
                double AIC = -2.0*ll + 2.0*num_params;
                double BIC = -2.0*ll + num_params*std::log((double)N);

                double criVal = (criterion_type==0 ? BIC : AIC);
                if(criVal < local_best){
                    local_best = criVal;
                    localModel.pos_start = pos_start;
                    localModel.pos_end   = pos_end;
                    localModel.cs        = cs;
                    localModel.ce        = ce;
                    localModel.beta[0]   = b0;
                    localModel.beta[1]   = b1;
                    localModel.logLik    = ll;
                    localModel.AIC       = AIC;
                    localModel.BIC       = BIC;
                }
            }
        }

        #pragma omp critical
        {
            if(local_best < bestCrit){
                bestCrit = local_best;
                bestM = localModel;
            }
        }
    } // end parallel

    return bestM;
}

// [[Rcpp::export]]
List runClusterSurvivalCovariates(IntegerVector deaths_r,
                                  IntegerVector pop_r,
                                  NumericMatrix X,
                                  int criterion_type=0)
{
    vector<long> deaths(deaths_r.size()), pop(pop_r.size());
    for(int i=0; i<deaths_r.size(); i++){
        deaths[i] = deaths_r[i];
        pop[i]    = pop_r[i];
    }

    long N = deaths_r.size();
    int K = X.ncol();

    ModelParams best = find_best_model_with_covariates(deaths, pop, X,
                                                       0, N-1,
                                                       K,
                                                       criterion_type);

    // 返回
    List out;
    out["pos_start"] = best.pos_start + 1;
    out["pos_end"]   = best.pos_end   + 1;
    out["cs"]        = best.cs        + 1;
    out["ce"]        = best.ce        + 1;
    out["logLik"]    = best.logLik;
    out["AIC"]       = best.AIC;
    out["BIC"]       = best.BIC;

    NumericMatrix Beta(2, K+1);
    for(int g=0; g<2; g++){
        for(int j=0; j<(K+1); j++){
            Beta(g,j) = best.beta[g][j];
        }
    }
    out["beta"] = Beta;

    return out;
}
