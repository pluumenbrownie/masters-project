//! https://cxx.rs/tutorial.html#calling-a-c-function-from-rust

fn main() {
    cxx_build::bridge("tests/compatability_tests.rs")
        .file("cpp_functions/Complexity.cpp")
        .std("c++11")
        .opt_level(3)
        .compile("complexity");

    println!("cargo:rerun-if-changed=cpp_functions/Complexity.cpp");
}
