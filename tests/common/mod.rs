//! Common functions for integration tests

use std::path::PathBuf;

pub fn get_resource_dir() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("resources/");
    path
}

pub fn get_test_file(file: &str) -> PathBuf {
    let mut path = get_resource_dir();
    path.push(file);
    path
}
