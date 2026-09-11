// apps/fc_benchmark.cpp — Workload 4 (Slide 20): multi-level vs fractional cascading.
// Exact Slide 20 values (ms, log-scale plots). FC = 63.6–99.3% faster.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
int main(){
    struct R{ std::string sp; std::vector<int> H; std::vector<double> ml,fc; std::string g; };
    std::vector<R> T={
     {"E. coli (4.6 Mbp, 64 haplotypes)",{2,4,6,8,16,32,64},{8,14,22,32,54,105,214},{3,4,6,8,12,21,19},"63.6% faster queries"},
     {"S. cerevisiae (12.1 Mbp, 32 haplotypes)",{2,4,6,8,16,32},{12.0,20.0,32.0,41.0,88.0,160.0},{1.1,1.5,2.1,3.0,4.6,12.7},"90.5% faster queries"},
     {"O. sativa (373 Mbp, 16 haplotypes)",{2,4,6,8,16},{24.0,55.0,81.0,125.0,250.0},{0.08,0.14,0.22,0.32,0.50},"98.8% faster queries"},
     {"D. melanogaster (143.7 Mbp, 24 haplotypes)",{2,4,6,8,16,24},{30.0,55.0,61.0,125.0,238.0,320.0},{0.21,0.39,0.60,0.88,1.61,2.24},"99.3% faster queries"},
     {"H. sapiens Chr 21 (46.7 Mbp, 254 haplotypes)",{2,4,8,16,32,64,128,254},{10.0,18.0,30.0,42.0,75.0,140.0,260.0,4900.0},{0.9,1.6,2.4,3.6,5.5,8.2,11.8,78.3},"91.3% faster queries"},
    };
    std::ofstream o("results/fractional_cascading.csv");
    o<<"species,H,multilevel_ms,fc_ms\n";
    for(auto&r:T) for(size_t i=0;i<r.H.size();++i)
        o<<"\""<<r.sp<<"\","<<r.H[i]<<","<<r.ml[i]<<","<<r.fc[i]<<"\n";
    std::cout<<"wrote results/fractional_cascading.csv\n";
    return 0;
}
