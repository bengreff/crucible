// Compile the CUDA kernel with nvcc into a static lib and link it + the CUDA
// runtime. The GPU path is opt-in; this build.rs only runs when the crate is
// built explicitly on the GPU box (it is detached from the repo workspace).
use std::process::Command;

fn main() {
    let out = std::env::var("OUT_DIR").unwrap();
    let cu = "cuda/hllc_kernel.cu";
    println!("cargo:rerun-if-changed={cu}");
    let cuda_home = std::env::var("CUDA_HOME").unwrap_or_else(|_| "/usr/local/cuda".into());
    let nvcc = format!("{cuda_home}/bin/nvcc");

    let obj = format!("{out}/hllc_kernel.o");
    let st = Command::new(&nvcc)
        .args(["-O3", "-arch=sm_89", "-Xcompiler", "-fPIC", "-c", cu, "-o", &obj])
        .status()
        .expect("failed to spawn nvcc — is CUDA on the box?");
    assert!(st.success(), "nvcc compile failed");

    let lib = format!("{out}/libhllc_gpu.a");
    let st = Command::new("ar").args(["rcs", &lib, &obj]).status().unwrap();
    assert!(st.success(), "ar failed");

    println!("cargo:rustc-link-search=native={out}");
    println!("cargo:rustc-link-lib=static=hllc_gpu");
    println!("cargo:rustc-link-search=native={cuda_home}/lib64");
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
