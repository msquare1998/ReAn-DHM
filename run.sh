# ============================================================================================
#  Re-An method for measuring lnZ of the 1D dimerzied Heisenberg model
#  Reference: https://arxiv.org/abs/2403.08642
#  @Yiming_Ding, Westlake University
#  Last updated: Mar 19, 2024
# ============================================================================================
L=6             # even number required in this program
epsilon=1e-1    # the magnitude of each ratio required
BETA=20         # inverse temperature
N_THM=5000      # thermalization time
N_STAT=10000    # measurement time
J_W0=1e-8       # minimum value of Jw0
J_W=1           # start from Jw = 1
nThread=10      # number of threads allocated
nBins=3         # number of bins

# =======================
#  For making diff dirs
# =======================
l_f="L"
epsilon_f="_epslion"
jw_f="_Jw"
jw0_f="_JwRef"
beta_f="_beta"
Lambda=`awk 'BEGIN{printf "%.0f\n",(1 * '$L')}'`
dir="$l_f$L$epsilon_f$epsilon$beta_f$BETA$jw0_f$J_W0$jw_f$J_W"

# =======================
#  Compile and run
# =======================
rm -rf data && mkdir data
rm -rf build && mkdir build
cd build && mkdir "../data/$dir"
mkdir "$dir" && cd "$dir"
cmake ../../ && make
./run $L $epsilon $BETA $Lambda $N_THM $N_STAT $J_W0 $J_W $nThread $nBins