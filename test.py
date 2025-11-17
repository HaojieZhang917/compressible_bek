# pep_quadratic.py
# 需在已安装 petsc4py 和 slepc4py 的环境下运行 （可用 mpirun -n 1 python pep_quadratic.py）
from petsc4py import PETSc
from slepc4py import SLEPc
import numpy as np

def build_1d_laplacian(N):
    """返回 PETSc 矩阵 K（1D Laplacian，N x N）"""
    K = PETSc.Mat().createAIJ([N, N], nnz=3)
    K.setFromOptions()
    K.setUp()
    for i in range(N):
        if i > 0:
            K[i, i-1] = -1.0
        K[i, i] = 2.0
        if i < N-1:
            K[i, i+1] = -1.0
    K.assemblyBegin()
    K.assemblyEnd()
    return K

def build_tridiag(N, diag=0.1, off= -0.05):
    """返回 PETSc 三对角矩阵 C"""
    C = PETSc.Mat().createAIJ([N, N], nnz=3)
    C.setFromOptions()
    C.setUp()
    for i in range(N):
        if i > 0:
            C[i, i-1] = off
        C[i, i] = diag
        if i < N-1:
            C[i, i+1] = off
    C.assemblyBegin()
    C.assemblyEnd()
    return C

def build_identity(N):
    M = PETSc.Mat().createAIJ([N, N], nnz=1)
    M.setFromOptions()
    M.setUp()
    for i in range(N):
        M[i, i] = 1.0
    M.assemblyBegin()
    M.assemblyEnd()
    return M

def main():
    # 初始化（通常 import slepc4py 会自动初始化 PETSc/SLEPc）
    # 建立矩阵
    N = 40  # 自己改
    K = build_1d_laplacian(N)
    C = build_tridiag(N, diag=0.05, off=-0.02)
    M = build_identity(N)

    # PEP 求解器
    pep = SLEPc.PEP().create()
    # 注意：SLEPc 中 PEP 的矩阵顺序通常是 A0 + lambda*A1 + lambda^2*A2
    # 所以对于 (lambda^2 M + lambda C + K) x = 0, 应传 [K, C, M]
    mats = [K, C, M]
    pep.setOperators(mats)
    pep.setProblemType(SLEPc.PEP.ProblemType.GENERAL)  # 一般多项式问题
    # 请求的特征值个数
    pep.setDimensions(nev=6)  # 要求 6 个特征值
    # 选择求解策略，例如靠近原点（INTERIOR），或最大模等，可用命令行覆盖
    pep.setWhichEigenpairs(SLEPc.PEP.Which.LARGEST_MAGNITUDE)

    # 从命令行读取可选项（例如 -pep_type TOAR 或 -pep_nev 10 等）
    pep.setFromOptions()

    # 求解
    pep.solve()

    nconv = pep.getConverged()
    print(f"Converged eigenpairs: {nconv}")

    vr = PETSc.Vec().createSeq(N)
    vi = PETSc.Vec().createSeq(N)
    for i in range(min(nconv, 6)):
        k = pep.getEigenpair(i, vr, vi)  # 返回复特征值（PetscScalar），向量放到 vr/vi
        # k 可能是复数（如果构建为复数算子或求到复特征值）
        try:
            val = complex(k)
        except Exception:
            # 有些版本返回 PetscScalar 对象，直接尝试 real/imag
            val = k
        print(f"#{i}: eigenvalue = {val.real} + {val.imag}j")
        # 可计算残差 ||P(lambda)x|| 来检验
        # 计算 y = A0 x + lambda A1 x + lambda^2 A2 x
        x = vr.copy()
        if vi.norm() > 1e-16:
            # 如果有虚部则 x 是复向量；这里仅打印实部范数作为示例
            pass
        # 计算残差（实值近似）
        tmp = PETSc.Vec().createSeq(N)
        # tmp = K*x
        K.mult(vr, tmp)
        # tmp += lambda*C*x
        tp2 = PETSc.Vec().createSeq(N)
        C.mult(vr, tp2)
        tmp.axpy(val, tp2)  # tmp += val * (C*x)
        # add lambda^2 * M * x
        tp3 = PETSc.Vec().createSeq(N)
        M.mult(vr, tp3)
        tmp.axpy(val*val, tp3)
        res_norm = tmp.norm()
        print(f"    residual norm (approx) = {res_norm:.3e}")

if __name__ == "__main__":
    main()
