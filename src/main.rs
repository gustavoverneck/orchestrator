use std::{fs, vec};
use std::process::Command;
use std::path::Path;
use rayon::prelude::*;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use console::style;

const TEMPLATE_FILE: &str = "input/template.f";
const TOV_FILE: &str = "input/tov.f";
const EOS_GEN_EXE: &str = "eos.exe";
const TOV_SOLVER_EXE: &str = "tov.exe";

#[derive(Debug)]
struct SimulationParam {
    model: String,
    b: String,
    x: String,
}

fn compile_fortran_codes() -> std::io::Result<()> {
    // Compile template.f -> eos_gen
    println!("Compiling {} to {}...", TEMPLATE_FILE, EOS_GEN_EXE);
    let status_eos = Command::new("gfortran")
        .arg(TEMPLATE_FILE)
        .arg("-o")
        .arg(EOS_GEN_EXE)
        .status()?;

    if !status_eos.success() {
        return Err(std::io::Error::new(std::io::ErrorKind::Other, "Failed to compile template.f"));
    }

    // Compile tov.f -> tov_solver
    println!("Compiling {} to {}...", TOV_FILE, TOV_SOLVER_EXE);
    let status_tov = Command::new("gfortran")
        .arg(TOV_FILE)
        .arg("-o")
        .arg(TOV_SOLVER_EXE)
        .status()?;

    if !status_tov.success() {
        return Err(std::io::Error::new(std::io::ErrorKind::Other, "Failed to compile tov.f"));
    }

    Ok(())
}

fn run_simulation(param: &SimulationParam, root_dir: &Path, pb: &ProgressBar) -> std::io::Result<()> {
    pb.set_message(format!("{} {:?}...", style("Init").cyan().bold(), param.x));
    let dir_path = Path::new("output").join(&param.model).join(&param.b).join(&param.x);
    fs::create_dir_all(&dir_path)?;

    // Paths to executables (must be absolute or relative to the new dir)
    // Using absolute paths is safer.
    let eos_gen_path = root_dir.join(EOS_GEN_EXE);
    let tov_solver_path = root_dir.join(TOV_SOLVER_EXE);
    // let plot_script_path = root_dir.join(PLOT_SCRIPT);

    // 1. Run EOS Generator
    pb.set_message(format!("{} {:?}...", style("EOS Gen").yellow().bold(), param.x));
    let run_eos = Command::new(&eos_gen_path)
        .arg(&param.x)
        .arg(&param.b)
        .arg(&param.model)
        .current_dir(&dir_path)
        .output()?; // Use output to suppress stdout unless error, or status to show it.

    if !run_eos.status.success() {
        pb.set_message(format!("{}", style("EOS Failed").red().bold()));
        eprintln!("EOS Generation failed in {:?}: {:?}", dir_path, String::from_utf8_lossy(&run_eos.stderr));
        return Ok(());
    }
    
    // Append terminator to eos.dat
    let eos_file = dir_path.join("eos.dat");
    if eos_file.exists() {
        use std::io::Write;
        let mut file = fs::OpenOptions::new()
            .write(true)
            .append(true)
            .open(eos_file)?;
        writeln!(file, "-1. -1. -1.")?;
    } else {
        eprintln!("Warning: eos.dat not found in {:?}", dir_path);
    }

    // 2. Run TOV Solver
    pb.set_message(format!("{} {:?}...", style("TOV Solver").blue().bold(), param.x));
    let run_tov = Command::new(&tov_solver_path)
        .current_dir(&dir_path)
        .output()?;

    if !run_tov.status.success() {
        pb.set_message(format!("{}", style("TOV Failed").red().bold()));
        eprintln!("TOV Solver failed in {:?}: {:?}", dir_path, String::from_utf8_lossy(&run_tov.stderr));
    }

    // 3. Run Plot Script
    // if plot_script_path.exists() {
    //     pb.set_message(format!("{} {:?}...", style("Plotting").magenta().bold(), param.x));
        
    //     let run_plot = Command::new("python")
    //         .arg(&plot_script_path) // Pass absolute path to script
    //         .current_dir(&dir_path)
    //         .output()?;
            
    //     if !run_plot.status.success() {
    //          pb.set_message(format!("{}", style("Plot Failed").red().bold()));
    //     }
    // }
    
    pb.inc(1); // Increment main progress if we tracked total, but here we update text
    Ok(())
}

fn get_x_params(
    negative_exp: bool,
    add_zero: bool,
    ei: i32,
    ef: i32,
    start: f64,
    end: f64,
    step: f64,
) -> Vec<String> {
    let mut params = Vec::new();

    if add_zero {
        params.push("0.0d0".to_string());
    }

    let epsilon = 1e-10;
    
    // But let's stick to the requested logic: step iteration.
    let mut j = start;
    while j < end - epsilon {
        let mut j_str = format!("{:.10}", j); // Avoid super long precisions like 0.30000000000000004
        // Trim trailing zeros but keep .0
        j_str = j_str.trim_end_matches('0').to_string();
        if j_str.ends_with('.') {
            j_str.push('0');
        }
        
        // If negative_exp is true
        if negative_exp {
            for e in ei..ef {
                params.push(format!("{}d-{}", j_str, e));
            }
        } else {
            for e in ei..ef {
                params.push(format!("{}d{}", j_str, e));
            }
        }
        j += step;
    }

    params
}

fn collect_results() -> std::io::Result<()> {
    let results_dir = Path::new("results");
    if !results_dir.exists() {
        fs::create_dir_all(results_dir)?;
    }
    let output_file = results_dir.join("results_aggregate.csv");
    let mut wtr = std::fs::File::create(&output_file)?;
    // Header
    use std::io::Write;
    writeln!(wtr, "model,b,csi,mmax,rmax")?;

    let output_dir = Path::new("output");
    if !output_dir.exists() {
        return Ok(());
    }

    // output/MODEL/B/CSI/tov.out
    for model_entry in fs::read_dir(output_dir)? {
        let model_entry = model_entry?;
        if !model_entry.path().is_dir() { continue; }
        let model_name = model_entry.file_name().into_string().unwrap();

        for b_entry in fs::read_dir(model_entry.path())? {
            let b_entry = b_entry?;
            if !b_entry.path().is_dir() { continue; }
            let b_name = b_entry.file_name().into_string().unwrap();

            for csi_entry in fs::read_dir(b_entry.path())? {
                let csi_entry = csi_entry?;
                if !csi_entry.path().is_dir() { continue; }
                let csi_name = csi_entry.file_name().into_string().unwrap();

                let tov_path = csi_entry.path().join("tov.out");
                if tov_path.exists() {
                    // Manually parse tov.out to find max mass
                    // tov.out columns: e0, m, r, mb
                    let content = fs::read_to_string(&tov_path)?;
                    let mut max_m = -1.0;
                    let mut r_at_max = -1.0;

                    for line in content.lines() {
                        let parts: Vec<&str> = line.split_whitespace().collect();
                        if parts.len() >= 3 {
                            if let (Ok(m), Ok(r)) = (parts[1].parse::<f64>(), parts[2].parse::<f64>()) {
                                if m > max_m {
                                    max_m = m;
                                    r_at_max = r;
                                }
                            }
                        }
                    }

                    if max_m > 0.0 {
                        writeln!(wtr, "{},{},{},{},{}", model_name, b_name, csi_name, max_m, r_at_max)?;
                    }
                }
            }
        }
    }
    println!("Aggregated results written to {:?}", output_file);
    Ok(())
}

fn main() -> std::io::Result<()> {
    // Check command line arguments for modes
    let args: Vec<String> = std::env::args().collect();
    let collect_only = args.contains(&"--collect-only".to_string());

    if !collect_only {
        // Defining parameter space
        let mut x_params = Vec::new();

        // Full range enabled
        x_params.extend(get_x_params(true, false, 12, 16, 1.0, 10.0, 1.0));
        x_params.extend(get_x_params(true, false, 9, 12, 1.0, 10.0, 0.1));
        x_params.push("0.0d0".to_string());

        let model_params = vec!["GM1", "GM3"];
        let b_params = vec!["1.0d15", "1.0d16", "1.0d17"];

        // Compile once
        if let Err(e) = compile_fortran_codes() {
            eprintln!("Detailed error compiling Fortran codes: {}", e);
            return Err(e);
        }
        
        let root_dir = std::env::current_dir()?;
        println!("Root dir: {:?}", root_dir);

        let mut params: Vec<SimulationParam> = Vec::new();

        for model in &model_params {
            for b in &b_params {
                for x in &x_params {
                    let param = SimulationParam {
                        model: model.to_string(),
                        b: b.to_string(),
                        x: x.clone(),
                    };
                    params.push(param);
                }
            }
        }

        println!("Total simulations generated: {}", params.len());

        // Setup worker bars for each thread
        let num_threads = rayon::current_num_threads();
        println!("Running on {} threads", num_threads);

        let multi_progress = MultiProgress::new();
        let main_pb = multi_progress.add(ProgressBar::new(params.len() as u64));
        main_pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"));
        main_pb.set_message("Total Progress");
        
        let worker_pbs: Vec<ProgressBar> = (0..num_threads).map(|i| {
            let pb = multi_progress.add(ProgressBar::new(0)); // Spinner
            pb.set_style(ProgressStyle::with_template(&format!("{{spinner}} [Thread {}] {{msg}}", i)).unwrap());
            pb
        }).collect();

        // Use Rayon for parallel processing
        params.par_iter().for_each(|param| {
            // Identify which thread we are on
            let thread_idx = rayon::current_thread_index().unwrap_or(0);
            // Select the corresponding progress bar
            let pb = &worker_pbs[thread_idx % worker_pbs.len()];

            if let Err(_e) = run_simulation(param, &root_dir, pb) {
                pb.set_message(format!("{} {:?}", style("Error").red().bold(), param.x));
            } else {
                // We don't "finish" these bars, just update message
                pb.set_message(format!("{} {:?}", style("Done").green().bold(), param.x));
            }
            main_pb.inc(1);
        });
        
        // Cleanup worker bars
        for pb in worker_pbs {
            pb.finish_and_clear();
        }
        
        main_pb.finish_with_message("All simulations completed.");
    } else {
        println!("Skipping simulations (--collect-only detected).");
    }

    // Aggregation Step
    println!("Collecting results...");
    if let Err(e) = collect_results() {
        eprintln!("Failed to collect results: {}", e);
    }

    // Run Plotting Script
    println!("Running plotting script with uv...");
    
    // Try 'uv run' which handles dependencies automatically via the script header
    let status_plot = Command::new("uv")
        .arg("run")
        .arg("input/plotting.py")
        .status();

    let msg = match status_plot {
        Ok(status) => {
             if status.success() { "Plotting successful." } else { "Plotting script returned error." }
        },
        Err(_) => "Failed to execute 'uv'. Make sure uv is installed and in PATH.",
    };
    println!("{}", msg);

    Ok(())
}


