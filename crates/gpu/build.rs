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
        "cuda/residency_class_r.cu",
        "cuda/residency_geometry.cu",
        "cuda/residency_engine.cu",
    ];
    let mut objs = Vec::new();
    for cu in sources {
        println!("cargo:rerun-if-changed={cu}");
        let stem = std::path::Path::new(cu)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap();
        let obj = format!("{out}/{stem}.o");
        // CRUCIBLE_NVCC_FMAD=0 builds the kernels without FMA contraction
        // (IEEE per-operation arithmetic, as the CPU reference): the
        // cross-device-ECT experiment (S14). Default: nvcc's contraction on.
        let fmad = std::env::var("CRUCIBLE_NVCC_FMAD").map_or(true, |v| v != "0");
        println!("cargo:rerun-if-env-changed=CRUCIBLE_NVCC_FMAD");
        let mut args = vec!["-O3", "-arch=sm_89", "-Xcompiler", "-fPIC"];
        if !fmad {
            args.push("-fmad=false");
        }
        args.extend_from_slice(&["-c", cu, "-o", &obj]);
        let st = Command::new(&nvcc)
            .args(&args)
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
