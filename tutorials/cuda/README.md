Incremental improvements on CUDA reduction algorithm.

Resources:
* CUDA programming guide: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#
* Old reduction example (some things are outdated): https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
    - NOTE: As of Volta (compute_70) instructions within a warp are no longer guaranteed to be synchronous, so `__syncwarp()` or the `shfl` instructions are needed to unroll the last warp (slide 21)
* Shared memory: https://developer.nvidia.com/blog/using-shared-memory-cuda-cc/
* Vectorised memory loads: https://developer.nvidia.com/blog/cuda-pro-tip-increase-performance-with-vectorized-memory-access/
* CUDA 9 warp level primitives (Volta): https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/
* Shuffle instructions: https://developer.nvidia.com/blog/faster-parallel-reductions-kepler/
* Cooperative groups: https://developer.nvidia.com/blog/cooperative-groups/
