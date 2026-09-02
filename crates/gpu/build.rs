// Compile the CUDA kernels with nvcc into a static lib and link it + the CUDA
// runtime. The GPU path is opt-in; this build.rs only runs when the crate is
// built explicitly on the GPU box (it is detached from the repo workspace).
use std::process::Command;

fn main() {
    let out = std::env::var("OUT_DIR").unwrap();
    let cuda_home = std::env::var("CUDA_HOME").unwrap_or_else(|_| "/usr/local/cuda".into());
    let nvcc = format!("{cuda_home}/bin/nvcc");

    // Each .cu -> its own object; all archived into one static lib.
    let sources = [
        "cuda/hllc_kernel.cu",
        "cuda/residency.cu",
        "cuda/residency_diffusion.cu",
        "cuda/residency_eos.cu",
        "cuda/residency_combustion.cu",
        "cuda/residency_blend_eos.cu",
    ];
    let mut objs = Vec::new();
    for cu in sources {
        println!("cargo:rerun-if-changed={cu}");
        let stem = std::path::Path::new(cu).file_stem().unwrap().to_str().unwrap();
        let obj = format!("{out}/{stem}.o");
        let st = Command::new(&nvcc)
            .args(["-O3", "-arch=sm_89", "-Xcompiler", "-fPIC", "-c", cu, "-o", &obj])
            .status()
            .expect("failed to spawn nvcc — is CUDA on the box?");
        assert!(st.success(), "nvcc compile failed for {cu}");
        objs.push(obj);
    }

    let lib = format!("{out}/libcrucible_gpu.a");
    let mut ar = Command::new("ar");
    ar.args(["rcs", &lib]);
    for o in &objs {
        ar.arg(o);
    }
    assert!(ar.status().unwrap().success(), "ar failed");

    println!("cargo:rustc-link-search=native={out}");
    println!("cargo:rustc-link-lib=static=crucible_gpu");
    println!("cargo:rustc-link-search=native={cuda_home}/lib64");
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
