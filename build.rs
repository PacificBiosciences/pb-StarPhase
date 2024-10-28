use std::error::Error;
use vergen::EmitBuilder;

fn main() -> Result<(), Box<dyn Error>> {
    EmitBuilder::builder()
        .fail_on_error()
        .all_git()
        .git_describe(true, false, Some("ThisPatternShouldNotMatchAnythingEver"))
        .emit()?;

    // emit build handles the git configuration and build.rs, but we also need to track the toml and src folder 
    let rerun_if_changed = "cargo:rerun-if-changed=Cargo.toml
cargo:rerun-if-changed=src";
    println!("{rerun_if_changed}");

    Ok(())
}