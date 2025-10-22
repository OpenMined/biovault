use crate::cli::examples;
use crate::error::Result;
use crate::project_spec::{self, InputSpec, OutputSpec, ParameterSpec, ProjectSpec};
use crate::types::InboxSubmission;
use colored::Colorize;
use dialoguer::{theme::ColorfulTheme, Confirm, Input, MultiSelect, Select};
use serde_yaml::Value;
use std::collections::HashSet;
use std::fs;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::process::Command;
use walkdir::WalkDir;

trait DialoguerResultExt<T> {
    fn cli_result(self) -> Result<T>;
}

impl<T> DialoguerResultExt<T> for std::result::Result<T, dialoguer::Error> {
    fn cli_result(self) -> Result<T> {
        self.map_err(|err| crate::error::Error::Anyhow(anyhow::anyhow!(err)))
    }
}

pub fn list_examples() -> Result<()> {
    let available = examples::get_available_examples()?;

    if available.is_empty() {
        println!("No example templates available.");
        return Ok(());
    }

    println!("Available example templates:\n");
    for (key, info) in available.iter() {
        println!("  {} - {}", key.green(), info.description);
        println!("    Template: {}", info.template);
        println!("    Files: {} files", info.files.len());
        println!();
    }

    println!("To create a project from an example:");
    println!("  bv project create --example <example-name>");

    Ok(())
}

pub fn view(path: Option<String>) -> Result<()> {
    let project_dir = path.unwrap_or_else(|| ".".to_string());
    let project_path = Path::new(&project_dir);

    if !project_path.exists() {
        return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
            "Project directory does not exist: {}",
            project_dir
        )));
    }

    let spec_path = project_path.join("project.yaml");
    if !spec_path.exists() {
        return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
            "project.yaml not found in {}",
            project_dir
        )));
    }

    let spec = project_spec::ProjectSpec::load(&spec_path)?;
    print_project_view(&spec, &project_dir);

    Ok(())
}

pub async fn create(
    name: Option<String>,
    folder: Option<String>,
    example: Option<String>,
    spec: Option<String>,
    input_to: Option<String>,
    output_from: Option<String>,
) -> Result<()> {
    // Determine which example to use
    let selected_example = if let Some(ex) = example {
        // Validate the example exists
        let available = examples::list_examples();
        if !available.contains(&ex) {
            return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
                "Unknown example '{}'. Available examples: {}",
                ex,
                available.join(", ")
            )));
        }
        Some(ex)
    } else {
        None
    };

    // Determine project name
    let mut prompted_for_name = false;
    let project_name = if let Some(ref ex) = selected_example {
        ex.clone()
    } else if let Some(n) = name {
        n
    } else if std::io::stdin().is_terminal() {
        prompted_for_name = true;
        Input::with_theme(&ColorfulTheme::default())
            .with_prompt("Project name")
            .validate_with(|input: &String| {
                if input.trim().is_empty() {
                    Err("Project name cannot be empty")
                } else {
                    Ok(())
                }
            })
            .interact_text()
            .cli_result()?
            .trim()
            .to_string()
    } else {
        println!("Enter project name:");
        let mut input = String::new();
        std::io::stdin().read_line(&mut input)?;
        input.trim().to_string()
    };

    // Determine folder path (default to ./<name>)
    let project_folder = folder.unwrap_or_else(|| format!("./{}", project_name));
    let project_path = Path::new(&project_folder);

    if project_path.exists() {
        return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
            "Project folder already exists: {}",
            project_folder
        )));
    }

    // Create base folder
    fs::create_dir_all(project_path)?;

    if let Some(ref example_name) = selected_example {
        // Use the new examples system to write files
        examples::write_example_to_directory(example_name, project_path)
            .map_err(crate::error::Error::Anyhow)?;

        // Load the generated spec to display summary
        let spec_path = project_path.join("project.yaml");
        let spec = ProjectSpec::load(&spec_path)?;
        println!("   (from example '{}')", example_name);
        print_project_summary(&spec, &project_folder);
    } else {
        // Load user email from config (if available)
        let email = match crate::config::Config::load() {
            Ok(cfg) => cfg.email,
            Err(_) => std::env::var("SYFTBOX_EMAIL").unwrap_or_else(|_| "".to_string()),
        };

        let author_default = if email.trim().is_empty() {
            "user@example.com".to_string()
        } else {
            email.trim().to_string()
        };

        // Only use interactive wizard if we have a terminal AND not in CI/non-interactive mode
        let terminal_available = std::io::stdin().is_terminal()
            && std::io::stdout().is_terminal()
            && std::env::var("CI").is_err()
            && std::env::var("BIOVAULT_NON_INTERACTIVE").is_err();

        // Load prepopulated inputs/outputs for pipeline composition
        // --output_from: new project's inputs FROM other project's outputs
        let prepopulated_inputs = if let Some(ref output_from_path) = output_from {
            let output_from_spec_path = Path::new(output_from_path).join("project.yaml");
            if !output_from_spec_path.exists() {
                return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
                    "Cannot load --output_from: project.yaml not found in {}",
                    output_from_path
                )));
            }
            let source_spec = project_spec::ProjectSpec::load(&output_from_spec_path)?;
            Some(source_spec.outputs)
        } else {
            None
        };

        // --input_to: new project's outputs TO other project's inputs
        let prepopulated_outputs = if let Some(ref input_to_path) = input_to {
            let input_to_spec_path = Path::new(input_to_path).join("project.yaml");
            if !input_to_spec_path.exists() {
                return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
                    "Cannot load --input_to: project.yaml not found in {}",
                    input_to_path
                )));
            }
            let source_spec = project_spec::ProjectSpec::load(&input_to_spec_path)?;
            Some(source_spec.inputs)
        } else {
            None
        };

        let spec_data = if let Some(spec_path) = spec {
            let spec_path = Path::new(&spec_path);
            let spec_loaded = project_spec::ProjectSpec::load(spec_path)?;
            project_spec::scaffold_from_spec(spec_loaded, project_path)?
        } else if terminal_available {
            let wizard_spec = run_project_spec_wizard(
                &project_name,
                &author_default,
                prompted_for_name,
                prepopulated_inputs,
                prepopulated_outputs,
            )?;
            project_spec::scaffold_from_spec(wizard_spec, project_path)?
        } else {
            // Non-interactive fallback - use prepopulated values if available
            let fallback_inputs = if let Some(ref prepop) = prepopulated_inputs {
                prepop
                    .iter()
                    .map(|output| InputSpec {
                        name: output.name.clone(),
                        raw_type: output.raw_type.clone(),
                        description: output.description.clone(),
                        format: output.format.clone(),
                        path: output.path.clone(),
                        mapping: None,
                    })
                    .collect()
            } else {
                Vec::new()
            };

            let fallback_outputs = if let Some(ref prepop) = prepopulated_outputs {
                prepop
                    .iter()
                    .map(|input| OutputSpec {
                        name: input.name.clone(),
                        raw_type: input.raw_type.clone(),
                        description: input.description.clone(),
                        format: input.format.clone(),
                        path: input.path.clone(),
                    })
                    .collect()
            } else {
                Vec::new()
            };

            let fallback = ProjectSpec {
                name: project_name.clone(),
                author: author_default.clone(),
                workflow: "workflow.nf".to_string(),
                template: Some("dynamic-nextflow".to_string()),
                version: Some("0.1.0".to_string()),
                assets: Vec::new(),
                parameters: Vec::new(),
                inputs: fallback_inputs,
                outputs: fallback_outputs,
            };
            project_spec::scaffold_from_spec(fallback, project_path)?
        };
        print_project_summary(&spec_data, &project_folder);
    }

    Ok(())
}

fn print_project_summary(spec: &ProjectSpec, folder: &str) {
    println!("\n✅ Created project '{}'", spec.name.bold());
    println!("   Location: {}\n", folder);

    // ASCII diagram
    let max_width = 60;

    // Helper to calculate visual width (emojis = 2 columns)
    let visual_width = |s: &str| -> usize {
        s.chars()
            .map(|c| if c as u32 > 0x1F300 { 2 } else { 1 })
            .sum()
    };

    println!("┌{}┐", "─".repeat(max_width - 2));
    println!("│{:^width$}│", spec.name, width = max_width - 2);
    println!("├{}┤", "─".repeat(max_width - 2));

    // Inputs
    if !spec.inputs.is_empty() {
        let count_str = spec.inputs.len().to_string();
        let header = format!("│ 📥 Inputs ({})", count_str);
        let padding = max_width - visual_width(&header) - 1;
        println!("{}{}│", header, " ".repeat(padding));
        for input in &spec.inputs {
            let line = format!("│   • {}: {}", input.name, input.raw_type);
            let padding = max_width - visual_width(&line) - 1;
            println!("{}{}│", line, " ".repeat(padding));
        }
    } else {
        let line = "│ 📥 Inputs: none";
        let padding = max_width - visual_width(line) - 1;
        println!("{}{}│", line, " ".repeat(padding));
    }

    println!("│{}│", " ".repeat(max_width - 2));

    // Outputs
    if !spec.outputs.is_empty() {
        let count_str = spec.outputs.len().to_string();
        let header = format!("│ 📤 Outputs ({})", count_str);
        let padding = max_width - visual_width(&header) - 1;
        println!("{}{}│", header, " ".repeat(padding));
        for output in &spec.outputs {
            let line = format!("│   • {}: {}", output.name, output.raw_type);
            let padding = max_width - visual_width(&line) - 1;
            println!("{}{}│", line, " ".repeat(padding));
        }
    } else {
        let line = "│ 📤 Outputs: none";
        let padding = max_width - visual_width(line) - 1;
        println!("{}{}│", line, " ".repeat(padding));
    }

    if !spec.parameters.is_empty() {
        println!("│{}│", " ".repeat(max_width - 2));
        let count_str = spec.parameters.len().to_string();
        let header = format!("│ ⚙️  Parameters ({})", count_str);
        let padding = max_width - visual_width(&header) - 1;
        println!("{}{}│", header, " ".repeat(padding));
        for param in &spec.parameters {
            let line = format!("│   • {}: {}", param.name, param.raw_type);
            let padding = max_width - visual_width(&line) - 1;
            println!("{}{}│", line, " ".repeat(padding));
        }
    }

    println!("└{}┘\n", "─".repeat(max_width - 2));

    // File TOC
    println!("📁 Files:");
    println!("   ├─ project.yaml");
    println!("   ├─ workflow.nf");
    if spec.assets.is_empty() {
        println!("   └─ assets/ (empty)");
    } else {
        println!(
            "   └─ assets/ ({} item{})",
            spec.assets.len(),
            if spec.assets.len() == 1 { "" } else { "s" }
        );
    }

    println!("\n💡 Next steps:");
    println!("   1. cd {}", folder);
    println!("   2. Edit workflow.nf to implement your logic");
    println!("   3. Run with: bv run . --<input_name> <path>");
}

fn print_project_view(spec: &ProjectSpec, folder: &str) {
    println!("\n📦 Project: {}", spec.name.bold());
    println!("   Location: {}\n", folder);

    // ASCII diagram
    let max_width = 60;

    // Helper to calculate visual width (emojis = 2 columns)
    let visual_width = |s: &str| -> usize {
        s.chars()
            .map(|c| if c as u32 > 0x1F300 { 2 } else { 1 })
            .sum()
    };

    println!("┌{}┐", "─".repeat(max_width - 2));
    println!("│{:^width$}│", spec.name, width = max_width - 2);
    println!("├{}┤", "─".repeat(max_width - 2));

    // Inputs
    if !spec.inputs.is_empty() {
        let count_str = spec.inputs.len().to_string();
        let header = format!("│ 📥 Inputs ({})", count_str);
        let padding = max_width - visual_width(&header) - 1;
        println!("{}{}│", header, " ".repeat(padding));
        for input in &spec.inputs {
            let line = format!("│   • {}: {}", input.name, input.raw_type);
            let padding = max_width - visual_width(&line) - 1;
            println!("{}{}│", line, " ".repeat(padding));
        }
    } else {
        let line = "│ 📥 Inputs: none";
        let padding = max_width - visual_width(line) - 1;
        println!("{}{}│", line, " ".repeat(padding));
    }

    println!("│{}│", " ".repeat(max_width - 2));

    // Outputs
    if !spec.outputs.is_empty() {
        let count_str = spec.outputs.len().to_string();
        let header = format!("│ 📤 Outputs ({})", count_str);
        let padding = max_width - visual_width(&header) - 1;
        println!("{}{}│", header, " ".repeat(padding));
        for output in &spec.outputs {
            let line = format!("│   • {}: {}", output.name, output.raw_type);
            let padding = max_width - visual_width(&line) - 1;
            println!("{}{}│", line, " ".repeat(padding));
        }
    } else {
        let line = "│ 📤 Outputs: none";
        let padding = max_width - visual_width(line) - 1;
        println!("{}{}│", line, " ".repeat(padding));
    }

    if !spec.parameters.is_empty() {
        println!("│{}│", " ".repeat(max_width - 2));
        let count_str = spec.parameters.len().to_string();
        let header = format!("│ ⚙️  Parameters ({})", count_str);
        let padding = max_width - visual_width(&header) - 1;
        println!("{}{}│", header, " ".repeat(padding));
        for param in &spec.parameters {
            let line = format!("│   • {}: {}", param.name, param.raw_type);
            let padding = max_width - visual_width(&line) - 1;
            println!("{}{}│", line, " ".repeat(padding));
        }
    }

    println!("└{}┘\n", "─".repeat(max_width - 2));

    // File TOC
    println!("📁 Files:");
    println!("   ├─ project.yaml");
    println!("   ├─ workflow.nf");
    if spec.assets.is_empty() {
        println!("   └─ assets/ (empty)");
    } else {
        println!(
            "   └─ assets/ ({} item{})",
            spec.assets.len(),
            if spec.assets.len() == 1 { "" } else { "s" }
        );
    }
}

pub async fn list(show_all: bool, show_full: bool) -> Result<()> {
    let inbox_path = get_inbox_path()?;

    if !inbox_path.exists() {
        println!("📭 No submissions in inbox");
        return Ok(());
    }

    let mut submissions = Vec::new();

    for entry in WalkDir::new(&inbox_path)
        .min_depth(2)
        .max_depth(2)
        .follow_links(false)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("yaml"))
    {
        let path = entry.path();

        let sender = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        match InboxSubmission::from_file(&path.to_path_buf()) {
            Ok(submission) => {
                if !show_all && submission.status == "rejected" {
                    continue;
                }

                let filename = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");

                let date_str = extract_date_from_filename(filename);

                submissions.push((sender.to_string(), submission, path.to_path_buf(), date_str));
            }
            Err(e) => {
                eprintln!("Warning: Failed to load submission from {:?}: {}", path, e);
            }
        }
    }

    if submissions.is_empty() {
        if !show_all {
            println!("📭 No active submissions in inbox (use --all to show rejected)");
        } else {
            println!("📭 No submissions in inbox");
        }
        return Ok(());
    }

    submissions.sort_by(|a, b| b.3.cmp(&a.3));

    let filtered_text = if !show_all { " active" } else { "" };
    println!(
        "📬 {}{} submission(s) in inbox:\n",
        submissions.len(),
        filtered_text
    );

    if show_full {
        for (i, (sender, submission, path, _date)) in submissions.iter().enumerate() {
            display_full_submission(i + 1, sender, submission, path);
        }
    } else {
        display_concise_list(&submissions);
    }

    Ok(())
}

fn prompt_non_empty(
    prompt: &str,
    default: Option<&str>,
    seen: &mut HashSet<String>,
) -> Result<String> {
    let theme = ColorfulTheme::default();
    loop {
        let mut builder = Input::with_theme(&theme).with_prompt(prompt);
        if let Some(d) = default {
            builder = builder.default(d.to_string());
        }
        let value: String = builder.interact_text().cli_result()?;
        let trimmed = value.trim();
        if trimmed.is_empty() {
            println!("Value cannot be empty.");
            continue;
        }
        if !seen.insert(trimmed.to_string()) {
            println!("'{}' already used, choose another.", trimmed);
            continue;
        }
        return Ok(trimmed.to_string());
    }
}

fn prompt_optional_string(prompt: &str) -> Result<Option<String>> {
    let theme = ColorfulTheme::default();
    let input: String = Input::with_theme(&theme)
        .with_prompt(prompt)
        .allow_empty(true)
        .interact_text()
        .cli_result()?;
    let trimmed = input.trim();
    if trimmed.is_empty() {
        Ok(None)
    } else {
        Ok(Some(trimmed.to_string()))
    }
}

fn run_project_spec_wizard(
    project_name: &str,
    default_author: &str,
    skip_name_prompt: bool,
    prepopulated_inputs: Option<Vec<OutputSpec>>,
    prepopulated_outputs: Option<Vec<InputSpec>>,
) -> Result<ProjectSpec> {
    let theme = ColorfulTheme::default();

    println!("\n🧪 Project wizard — let's describe your Nextflow wrapper");

    let name: String = if skip_name_prompt {
        project_name.to_string()
    } else {
        Input::with_theme(&theme)
            .with_prompt("Project name (identifier)")
            .default(project_name.to_string())
            .interact_text()
            .cli_result()?
    };

    let author: String = Input::with_theme(&theme)
        .with_prompt("Author (email)")
        .default(default_author.to_string())
        .interact_text()
        .cli_result()?;

    let workflow: String = Input::with_theme(&theme)
        .with_prompt("Workflow script filename (filepath)")
        .default("workflow.nf".to_string())
        .interact_text()
        .cli_result()?;

    let version: Option<String> = {
        let value: String = Input::with_theme(&theme)
            .with_prompt("Project version (semver)")
            .default("0.1.0".to_string())
            .interact_text()
            .cli_result()?;
        let trimmed = value.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    };

    println!("\n📥 Inputs define data this workflow receives");
    println!("Supported types: String, Bool, File, Directory, ParticipantSheet, GenotypeRecord, List[...], Map[String, ...], Optional(?)");

    // Prepopulate inputs from --output_from (convert OutputSpec -> InputSpec)
    let mut inputs = if let Some(ref prepop) = prepopulated_inputs {
        println!(
            "✨ Prepopulated {} input(s) from source project's outputs",
            prepop.len()
        );
        prepop
            .iter()
            .map(|output| InputSpec {
                name: output.name.clone(),
                raw_type: output.raw_type.clone(),
                description: output.description.clone(),
                format: output.format.clone(),
                path: output.path.clone(),
                mapping: None, // Outputs don't have mapping
            })
            .collect()
    } else {
        Vec::new()
    };
    let mut seen_inputs: HashSet<String> = inputs.iter().map(|i| i.name.clone()).collect();
    while Confirm::with_theme(&theme)
        .with_prompt(if inputs.is_empty() {
            "Add an input?"
        } else {
            "Add another input?"
        })
        .default(inputs.is_empty())
        .interact()
        .cli_result()?
    {
        let name = prompt_non_empty(
            "Input identifier (key name, e.g. 'rows', 'files')",
            None,
            &mut seen_inputs,
        )?;
        let raw_type = loop {
            let ty: String = Input::with_theme(&theme)
                .with_prompt("Input type (e.g. List[GenotypeRecord])")
                .interact_text()
                .cli_result()?;
            match project_spec::validate_type_expr(&ty) {
                Ok(_) => break ty.trim().to_string(),
                Err(e) => {
                    println!("{}", format!("Invalid type: {}", e).red());
                }
            }
        };
        let description = prompt_optional_string("Description (optional)")?;
        let format = if raw_type.contains("ParticipantSheet") {
            prompt_optional_string("Format hint (optional, e.g. csv, tsv)")?
        } else {
            None
        };
        let path = prompt_optional_string("Default path/pattern (optional filepath)")?;

        inputs.push(InputSpec {
            name,
            raw_type,
            description,
            format,
            path,
            mapping: None,
        });
    }

    println!("\n📤 Outputs define data this workflow produces");

    // Prepopulate outputs from --input_to (convert InputSpec -> OutputSpec)
    let mut outputs = if let Some(ref prepop) = prepopulated_outputs {
        println!(
            "✨ Prepopulated {} output(s) from source project's inputs",
            prepop.len()
        );
        prepop
            .iter()
            .map(|input| OutputSpec {
                name: input.name.clone(),
                raw_type: input.raw_type.clone(),
                description: input.description.clone(),
                format: input.format.clone(),
                path: input.path.clone(),
            })
            .collect()
    } else {
        Vec::new()
    };
    let mut seen_outputs: HashSet<String> = outputs.iter().map(|o| o.name.clone()).collect();
    while Confirm::with_theme(&theme)
        .with_prompt(if outputs.is_empty() {
            "Add an output?"
        } else {
            "Add another output?"
        })
        .default(outputs.is_empty())
        .interact()
        .cli_result()?
    {
        let name = prompt_non_empty(
            "Output identifier (key name, e.g. 'scored_sheet')",
            None,
            &mut seen_outputs,
        )?;
        let raw_type = loop {
            let ty: String = Input::with_theme(&theme)
                .with_prompt("Output type (e.g. File, ParticipantSheet)")
                .interact_text()
                .cli_result()?;
            match project_spec::validate_type_expr(&ty) {
                Ok(_) => break ty.trim().to_string(),
                Err(e) => {
                    println!("{}", format!("Invalid type: {}", e).red());
                }
            }
        };
        let description = prompt_optional_string("Description (optional)")?;
        let default_path = format!("results/{}", name.replace(':', "_"));
        let path = {
            let value: String = Input::with_theme(&theme)
                .with_prompt("Output filepath (relative to results/)")
                .default(default_path)
                .interact_text()
                .cli_result()?;
            let trimmed = value.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        };
        let format = if raw_type.contains("ParticipantSheet") {
            prompt_optional_string("Format hint (optional, e.g. csv, tsv)")?
        } else {
            None
        };

        outputs.push(OutputSpec {
            name,
            raw_type,
            description,
            format,
            path,
        });
    }

    println!("\n⚙️  Parameters define runtime configuration toggles (optional)");
    println!("These are exposed as context.params.* in your workflow");

    let mut parameters = Vec::new();
    let mut seen_params = HashSet::new();
    if Confirm::with_theme(&theme)
        .with_prompt("Add parameters?")
        .default(false)
        .interact()
        .cli_result()?
    {
        while Confirm::with_theme(&theme)
            .with_prompt("Add a parameter?")
            .default(parameters.is_empty())
            .interact()
            .cli_result()?
        {
            let name = prompt_non_empty(
                "Parameter identifier (key name, e.g. 'use_cache')",
                None,
                &mut seen_params,
            )?;
            let param_type_options = vec!["String", "Bool", "Enum"];
            let choice = Select::with_theme(&theme)
                .with_prompt("Parameter type")
                .items(&param_type_options)
                .default(0)
                .interact()
                .cli_result()?;
            let raw_type = param_type_options[choice].to_string();

            let description = prompt_optional_string("Description (optional)")?;

            let mut default_value: Option<Value> = None;
            let mut choices_list: Option<Vec<String>> = None;

            match raw_type.as_str() {
                "String" => {
                    if let Some(def) = prompt_optional_string("Default value (optional)")? {
                        default_value = Some(Value::String(def));
                    }
                }
                "Bool" => {
                    let default_bool = Confirm::with_theme(&theme)
                        .with_prompt("Default to true?")
                        .default(false)
                        .interact()
                        .cli_result()?;
                    default_value = Some(Value::Bool(default_bool));
                }
                "Enum" => {
                    let choices_input: String = Input::with_theme(&theme)
                        .with_prompt("Comma separated choices (e.g. a,b,c)")
                        .interact_text()
                        .cli_result()?;
                    let mut parsed: Vec<String> = choices_input
                        .split(',')
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .collect();
                    parsed.sort();
                    parsed.dedup();
                    if parsed.is_empty() {
                        return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
                            "Enum parameter must have at least one choice"
                        )));
                    }
                    let default_choice =
                        prompt_optional_string("Default choice (leave blank to skip)")?;
                    if let Some(ref def) = default_choice {
                        if !parsed.contains(def) {
                            return Err(crate::error::Error::Anyhow(anyhow::anyhow!(
                                "Default '{}' not one of {:?}",
                                def,
                                parsed
                            )));
                        }
                        default_value = Some(Value::String(def.clone()));
                    }
                    choices_list = Some(parsed);
                }
                _ => {}
            }

            let advanced_flag = if Confirm::with_theme(&theme)
                .with_prompt("Mark as advanced (hide from simple wizards)?")
                .default(false)
                .interact()
                .cli_result()?
            {
                Some(true)
            } else {
                None
            };

            parameters.push(ParameterSpec {
                name,
                raw_type,
                description,
                default: default_value,
                choices: choices_list,
                advanced: advanced_flag,
            });
        }
    }

    println!("\n📁 Assets are static files bundled with your workflow (optional)");

    let mut assets = Vec::new();
    if Confirm::with_theme(&theme)
        .with_prompt("Add asset files?")
        .default(false)
        .interact()
        .cli_result()?
    {
        loop {
            let asset: String = Input::with_theme(&theme)
                .with_prompt("Asset filepath (relative to assets/, blank to finish)")
                .allow_empty(true)
                .interact_text()
                .cli_result()?;
            let trimmed = asset.trim();
            if trimmed.is_empty() {
                break;
            }
            assets.push(trimmed.to_string());
        }
    }

    Ok(ProjectSpec {
        name,
        author,
        workflow,
        template: Some("dynamic-nextflow".to_string()),
        version,
        assets,
        parameters,
        inputs,
        outputs,
    })
}

fn display_concise_list(submissions: &[(String, InboxSubmission, PathBuf, String)]) {
    println!(
        "{:<4} {:<12} {:<20} {:<25} {:<15} {:<30}",
        "#".bold(),
        "Date".bold(),
        "Project".bold(),
        "From".bold(),
        "Participants".bold(),
        "ID".bold()
    );
    println!("{}", "-".repeat(110));

    for (i, (sender, submission, path, date)) in submissions.iter().enumerate() {
        let status_icon = match submission.status.as_str() {
            "pending" => "⏳",
            "approved" => "✅",
            "rejected" => "❌",
            "reviewing" => "🔍",
            _ => "❓",
        };

        let participants_count = submission
            .participants
            .as_ref()
            .map(|p| p.len().to_string())
            .unwrap_or_else(|| "0".to_string());

        let id = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        let project_name = if submission.name.len() > 18 {
            format!("{}...", &submission.name[..15])
        } else {
            submission.name.clone()
        };

        let sender_display = if sender.len() > 23 {
            format!("{}...", &sender[..20])
        } else {
            sender.clone()
        };

        let id_display = if id.len() > 28 {
            format!("{}...", &id[..25])
        } else {
            id.to_string()
        };

        println!(
            "{:<4} {:<12} {} {:<18} {:<25} {:<15} {:<30}",
            format!("{}", i + 1),
            date,
            status_icon,
            project_name,
            sender_display.cyan(),
            participants_count,
            id_display.dimmed()
        );
    }
}

fn display_full_submission(index: usize, sender: &str, submission: &InboxSubmission, path: &Path) {
    let status_icon = match submission.status.as_str() {
        "pending" => "⏳",
        "approved" => "✅",
        "rejected" => "❌",
        "reviewing" => "🔍",
        _ => "❓",
    };

    println!("{}. {} {}", index, status_icon, submission.name.bold());
    println!("   From: {}", sender.cyan());
    println!("   Author: {}", submission.author);
    println!("   Status: {}", format_status(&submission.status));

    if let Some(datasites) = &submission.datasites {
        if !datasites.is_empty() {
            println!("   Datasites: {}", datasites.join(", "));
        }
    }

    if let Some(participants) = &submission.participants {
        if !participants.is_empty() {
            println!("   Participants: {} participant(s)", participants.len());
            if participants.len() <= 3 {
                for participant in participants {
                    println!("     - {}", participant.dimmed());
                }
            }
        }
    }

    println!("   Syft URL: {}", submission.syft_url.dimmed());

    if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
        println!("   ID: {}", filename.dimmed());
    }

    println!();
}

fn extract_date_from_filename(filename: &str) -> String {
    if filename.len() < 19 {
        return "unknown".to_string();
    }

    let parts: Vec<&str> = filename.split('-').collect();
    if parts.len() < 4 {
        return "unknown".to_string();
    }

    let len = parts.len();

    if let (Ok(year), Ok(month), Ok(day)) = (
        parts[len - 4].parse::<u32>(),
        parts[len - 3].parse::<u32>(),
        parts[len - 2].parse::<u32>(),
    ) {
        if (2020..=2100).contains(&year) && (1..=12).contains(&month) && (1..=31).contains(&day) {
            return format!("{:04}-{:02}-{:02}", year, month, day);
        }
    }

    "unknown".to_string()
}

fn get_inbox_path() -> Result<PathBuf> {
    Ok(crate::config::get_biovault_home()?.join("inbox"))
}

fn format_status(status: &str) -> String {
    match status {
        "pending" => status.yellow().to_string(),
        "approved" => status.green().to_string(),
        "rejected" => status.red().to_string(),
        "reviewing" => status.blue().to_string(),
        _ => status.to_string(),
    }
}

// tests moved to bottom of file

pub async fn show(reference: &str, show_all: bool) -> Result<()> {
    let submissions = load_submissions(show_all)?;

    if submissions.is_empty() {
        println!("📭 No submissions found");
        return Ok(());
    }

    if let Ok(index) = reference.parse::<usize>() {
        if index > 0 && index <= submissions.len() {
            show_submission_detail(&submissions[index - 1]);
            return Ok(());
        } else {
            println!(
                "❌ Invalid index: {}. Valid range: 1-{}",
                index,
                submissions.len()
            );
            return Ok(());
        }
    }

    let reference_lower = reference.to_lowercase();
    let matches: Vec<_> = submissions
        .iter()
        .filter(|(_, submission, path, _)| {
            if submission.name.to_lowercase().contains(&reference_lower) {
                return true;
            }
            if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                if filename.to_lowercase().contains(&reference_lower) {
                    return true;
                }
            }
            false
        })
        .collect();

    match matches.len() {
        0 => {
            println!("❌ No submission found matching: {}", reference);
            println!("💡 Try using: index number (1,2,3...), partial hash, or project name");
        }
        1 => {
            show_submission_detail(matches[0]);
        }
        _ => {
            println!("⚠️  Multiple submissions match '{}':", reference);
            for (i, (sender, submission, path, _)) in matches.iter().enumerate() {
                let id = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                println!(
                    "  {}. {} from {} ({})",
                    i + 1,
                    submission.name.bold(),
                    sender.cyan(),
                    &id[..id.len().min(8)].dimmed()
                );
            }
            println!("\n💡 Use a more specific reference or the full hash");
        }
    }

    Ok(())
}

pub async fn interactive(show_all: bool) -> Result<()> {
    let mut show_rejected = show_all;
    let mut current_index = 0;

    loop {
        let submissions = load_submissions(show_rejected)?;

        if submissions.is_empty() {
            println!("📭 No submissions found");
            if !show_rejected {
                let toggle = Confirm::with_theme(&ColorfulTheme::default())
                    .with_prompt("Show rejected submissions?")
                    .default(false)
                    .interact()
                    .map_err(|e| anyhow::anyhow!("Confirmation error: {}", e))?;

                if toggle {
                    show_rejected = true;
                    continue;
                }
            }
            return Ok(());
        }

        let items: Vec<String> = submissions
            .iter()
            .map(|(sender, submission, path, date)| {
                let status_icon = match submission.status.as_str() {
                    "pending" => "⏳",
                    "approved" => "✅",
                    "rejected" => "❌",
                    "reviewing" => "🔍",
                    _ => "❓",
                };

                let id = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                let short_id = if id.len() > 8 {
                    &id[id.len() - 8..]
                } else {
                    id
                };

                format!(
                    "{} {} - {} from {} [{}]",
                    status_icon, submission.name, date, sender, short_id
                )
            })
            .collect();

        let mut menu_items = items.clone();
        menu_items.push(format!(
            "🔄 Toggle view (currently: {})",
            if show_rejected { "all" } else { "active only" }
        ));
        menu_items.push("⚡ Batch actions".to_string());
        menu_items.push("📤 Exit".to_string());

        let selection = Select::with_theme(&ColorfulTheme::default())
            .with_prompt(
                "📬 Select submission or action (↑↓ to navigate, Enter to select, Esc to exit)",
            )
            .items(&menu_items)
            .default(current_index.min(menu_items.len() - 1))
            .interact_opt()
            .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;

        match selection {
            Some(index) if index < submissions.len() => {
                current_index = index;
                println!();

                let selected_submission = &submissions[index];
                show_submission_detail(selected_submission);

                println!("\n{}", "─".repeat(80));

                let actions = vec![
                    "📋 Back to list",
                    "❌ Reject submission",
                    "🔍 Mark for review",
                    "🧪 Test with mock data",
                    "📁 Show files",
                ];

                let action = Select::with_theme(&ColorfulTheme::default())
                    .with_prompt("Choose action")
                    .items(&actions)
                    .default(0)
                    .interact_opt()
                    .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;

                match action {
                    Some(0) => continue,
                    Some(1) => {
                        reject_submission(selected_submission)?;
                        println!("✅ Submission marked as rejected");
                    }
                    Some(2) => {
                        review_submission(selected_submission)?;
                        println!("✅ Submission marked for review");
                    }
                    Some(3) => {
                        test_submission(selected_submission).await?;
                    }
                    Some(4) => {
                        show_submission_files(selected_submission)?;
                    }
                    _ => continue,
                }
            }
            Some(index) if index == submissions.len() => {
                show_rejected = !show_rejected;
                current_index = 0;
            }
            Some(index) if index == submissions.len() + 1 => {
                batch_actions(&submissions, show_rejected)?;
                current_index = 0;
            }
            _ => {
                println!("👋 Exiting inbox");
                break;
            }
        }
    }

    Ok(())
}

fn batch_actions(
    submissions: &[(String, InboxSubmission, PathBuf, String)],
    _show_rejected: bool,
) -> Result<()> {
    let items: Vec<String> = submissions
        .iter()
        .enumerate()
        .map(|(i, (sender, submission, _path, date))| {
            let status_icon = match submission.status.as_str() {
                "pending" => "⏳",
                "approved" => "✅",
                "rejected" => "❌",
                "reviewing" => "🔍",
                _ => "❓",
            };
            format!(
                "{}. {} {} - {} from {}",
                i + 1,
                status_icon,
                submission.name,
                date,
                sender
            )
        })
        .collect();

    let selections = MultiSelect::with_theme(&ColorfulTheme::default())
        .with_prompt("Select submissions for batch action (Space to select, Enter to confirm)")
        .items(&items)
        .interact_opt()
        .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;

    if let Some(indices) = selections {
        if indices.is_empty() {
            println!("No submissions selected");
            return Ok(());
        }

        let actions = vec![
            "❌ Reject selected",
            "🔍 Mark selected for review",
            "🧪 Test selected with mock data",
            "Cancel",
        ];

        let action = Select::with_theme(&ColorfulTheme::default())
            .with_prompt("Choose batch action")
            .items(&actions)
            .default(0)
            .interact_opt()
            .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;

        match action {
            Some(0) => {
                for &idx in &indices {
                    reject_submission(&submissions[idx])?;
                }
                println!("✅ {} submission(s) marked as rejected", indices.len());
            }
            Some(1) => {
                for &idx in &indices {
                    review_submission(&submissions[idx])?;
                }
                println!("✅ {} submission(s) marked for review", indices.len());
            }
            Some(2) => {
                for &idx in &indices {
                    println!("\nTesting submission {}...", idx + 1);
                    drop(test_submission(&submissions[idx]));
                }
            }
            _ => {
                println!("Batch action cancelled");
            }
        }
    }

    Ok(())
}

fn load_submissions(show_all: bool) -> Result<Vec<(String, InboxSubmission, PathBuf, String)>> {
    let inbox_path = get_inbox_path()?;

    if !inbox_path.exists() {
        return Ok(Vec::new());
    }

    let mut submissions = Vec::new();

    for entry in WalkDir::new(&inbox_path)
        .min_depth(2)
        .max_depth(2)
        .follow_links(false)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("yaml"))
    {
        let path = entry.path();

        let sender = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        match InboxSubmission::from_file(&path.to_path_buf()) {
            Ok(submission) => {
                if !show_all && submission.status == "rejected" {
                    continue;
                }

                let filename = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");

                let date_str = extract_date_from_filename(filename);

                submissions.push((sender.to_string(), submission, path.to_path_buf(), date_str));
            }
            Err(e) => {
                eprintln!("Warning: Failed to load submission from {:?}: {}", path, e);
            }
        }
    }

    submissions.sort_by(|a, b| b.3.cmp(&a.3));

    Ok(submissions)
}

fn show_submission_detail(submission_data: &(String, InboxSubmission, PathBuf, String)) {
    let (sender, submission, path, date) = submission_data;

    let status_icon = match submission.status.as_str() {
        "pending" => "⏳",
        "approved" => "✅",
        "rejected" => "❌",
        "reviewing" => "🔍",
        _ => "❓",
    };

    println!("{}", "═".repeat(80));
    println!("{} {}", status_icon, submission.name.bold().cyan());
    println!("{}", "─".repeat(80));

    println!("📅 Date:        {}", date.bold());
    println!("📧 From:        {}", sender.cyan());
    println!("✍️  Author:      {}", submission.author);
    println!("📊 Status:      {}", format_status(&submission.status));

    if let Some(datasites) = &submission.datasites {
        if !datasites.is_empty() {
            println!("🏢 Datasites:   {}", datasites.join(", "));
        }
    }

    if let Some(participants) = &submission.participants {
        if !participants.is_empty() {
            println!(
                "👥 Participants: {} participant(s)",
                participants.len().to_string().bold()
            );
            for (i, participant) in participants.iter().enumerate() {
                if i < 5 {
                    println!("     {}. {}", i + 1, participant.dimmed());
                } else if i == 5 {
                    println!("     ... and {} more", participants.len() - 5);
                    break;
                }
            }
        }
    } else {
        println!("👥 Participants: None specified");
    }

    println!("🔗 Syft URL:    {}", submission.syft_url.dimmed());

    if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
        println!("🆔 ID:          {}", filename.dimmed());

        if filename.len() > 8 {
            let short_id = &filename[filename.len() - 8..];
            println!("🏷️  Short ID:    {}", short_id.yellow());
        }
    }

    println!("{}", "═".repeat(80));
}

fn show_submission_files(
    submission_data: &(String, InboxSubmission, PathBuf, String),
) -> Result<()> {
    // no interactive input here to avoid blocking tests
    let (_sender, _submission, path, _date) = submission_data;

    let submission_dir = path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Invalid submission path"))?;

    println!("\n📁 Files in submission directory:");
    println!("{}", "─".repeat(40));

    for entry in WalkDir::new(submission_dir)
        .max_depth(3)
        .follow_links(false)
        .into_iter()
        .filter_map(|e| e.ok())
    {
        if entry.file_type().is_file() {
            let relative_path = entry
                .path()
                .strip_prefix(submission_dir)
                .unwrap_or(entry.path());
            println!("  📄 {}", relative_path.display());
        }
    }

    println!("{}", "─".repeat(40));

    Ok(())
}

pub async fn reject(reference: Option<String>) -> Result<()> {
    if let Some(ref_str) = reference {
        let submissions = load_submissions(true)?;

        if let Ok(index) = ref_str.parse::<usize>() {
            if index > 0 && index <= submissions.len() {
                reject_submission(&submissions[index - 1])?;
                println!("✅ Submission marked as rejected");
                return Ok(());
            }
        }

        let reference_lower = ref_str.to_lowercase();
        let matches: Vec<_> = submissions
            .iter()
            .filter(|(_, submission, path, _)| {
                if submission.name.to_lowercase().contains(&reference_lower) {
                    return true;
                }
                if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                    if filename.to_lowercase().contains(&reference_lower) {
                        return true;
                    }
                }
                false
            })
            .collect();

        match matches.len() {
            0 => {
                println!("❌ No submission found matching: {}", ref_str);
            }
            1 => {
                reject_submission(matches[0])?;
                println!("✅ Submission marked as rejected");
            }
            _ => {
                println!(
                    "⚠️  Multiple submissions match '{}'. Please be more specific.",
                    ref_str
                );
            }
        }
    } else {
        println!("❌ Please specify a submission reference (index, hash, or name)");
    }

    Ok(())
}

pub async fn review(reference: Option<String>) -> Result<()> {
    if let Some(ref_str) = reference {
        let submissions = load_submissions(true)?;

        if let Ok(index) = ref_str.parse::<usize>() {
            if index > 0 && index <= submissions.len() {
                review_submission(&submissions[index - 1])?;
                println!("✅ Submission marked for review");
                return Ok(());
            }
        }

        let reference_lower = ref_str.to_lowercase();
        let matches: Vec<_> = submissions
            .iter()
            .filter(|(_, submission, path, _)| {
                if submission.name.to_lowercase().contains(&reference_lower) {
                    return true;
                }
                if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                    if filename.to_lowercase().contains(&reference_lower) {
                        return true;
                    }
                }
                false
            })
            .collect();

        match matches.len() {
            0 => {
                println!("❌ No submission found matching: {}", ref_str);
            }
            1 => {
                review_submission(matches[0])?;
                println!("✅ Submission marked for review");
            }
            _ => {
                println!(
                    "⚠️  Multiple submissions match '{}'. Please be more specific.",
                    ref_str
                );
            }
        }
    } else {
        println!("❌ Please specify a submission reference (index, hash, or name)");
    }

    Ok(())
}

pub async fn test(reference: Option<String>) -> Result<()> {
    if let Some(ref_str) = reference {
        let submissions = load_submissions(true)?;

        if let Ok(index) = ref_str.parse::<usize>() {
            if index > 0 && index <= submissions.len() {
                test_submission(&submissions[index - 1]).await?;
                return Ok(());
            }
        }

        let reference_lower = ref_str.to_lowercase();
        let matches: Vec<_> = submissions
            .iter()
            .filter(|(_, submission, path, _)| {
                if submission.name.to_lowercase().contains(&reference_lower) {
                    return true;
                }
                if let Some(filename) = path.file_stem().and_then(|s| s.to_str()) {
                    if filename.to_lowercase().contains(&reference_lower) {
                        return true;
                    }
                }
                false
            })
            .collect();

        match matches.len() {
            0 => {
                println!("❌ No submission found matching: {}", ref_str);
            }
            1 => {
                test_submission(matches[0]).await?;
            }
            _ => {
                println!(
                    "⚠️  Multiple submissions match '{}'. Please be more specific.",
                    ref_str
                );
            }
        }
    } else {
        println!("❌ Please specify a submission reference (index, hash, or name)");
    }

    Ok(())
}

fn reject_submission(submission_data: &(String, InboxSubmission, PathBuf, String)) -> Result<()> {
    let (_sender, mut submission, path, _date) = submission_data.clone();

    submission.status = "rejected".to_string();

    let yaml = serde_yaml::to_string(&submission)?;
    fs::write(&path, yaml)?;

    Ok(())
}

fn review_submission(submission_data: &(String, InboxSubmission, PathBuf, String)) -> Result<()> {
    let (_sender, mut submission, path, _date) = submission_data.clone();

    submission.status = "reviewing".to_string();

    let yaml = serde_yaml::to_string(&submission)?;
    fs::write(&path, yaml)?;

    Ok(())
}

async fn test_submission(
    submission_data: &(String, InboxSubmission, PathBuf, String),
) -> Result<()> {
    let (_sender, submission, path, _date) = submission_data;

    println!("🧪 Testing submission: {}", submission.name.bold());

    let submission_dir = path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Invalid submission path"))?;

    let workflow_file = submission_dir.join("workflow.nf");
    if !workflow_file.exists() {
        println!("❌ No workflow.nf file found in submission");
        return Ok(());
    }

    if let Some(participants) = &submission.participants {
        if participants.is_empty() {
            println!("❌ No participants specified in submission");
            return Ok(());
        }

        let participant = if participants.len() == 1 {
            &participants[0]
        } else {
            println!("Select participant to test with:");
            let selection = Select::with_theme(&ColorfulTheme::default())
                .items(participants)
                .default(0)
                .interact()
                .map_err(|e| anyhow::anyhow!("Selection error: {}", e))?;
            &participants[selection]
        };

        println!("🔬 Running test with participant: {}", participant.cyan());

        let mock_data_path = crate::config::get_biovault_home()?
            .join("data")
            .join("sample")
            .join(participant);

        if !mock_data_path.exists() {
            println!("⚠️  No mock data found for participant: {}", participant);
            println!("   Expected path: {}", mock_data_path.display());

            let download = Confirm::with_theme(&ColorfulTheme::default())
                .with_prompt("Download sample data for this participant?")
                .default(true)
                .interact()
                .map_err(|e| anyhow::anyhow!("Confirmation error: {}", e))?;

            if download {
                println!("📥 Downloading sample data for {}...", participant);
                Command::new("bv")
                    .args(["sample-data", "fetch", participant])
                    .status()
                    .map_err(|e| anyhow::anyhow!("Failed to download sample data: {}", e))?;
            } else {
                return Ok(());
            }
        }

        let output_dir = submission_dir.join("test_output");
        fs::create_dir_all(&output_dir)?;

        println!("🚀 Running nextflow workflow...");
        println!("   Input: {}", mock_data_path.display());
        println!("   Output: {}", output_dir.display());

        let status = Command::new("nextflow")
            .arg("run")
            .arg(&workflow_file)
            .arg("--input")
            .arg(&mock_data_path)
            .arg("--output")
            .arg(&output_dir)
            .current_dir(submission_dir)
            .status()
            .map_err(|e| anyhow::anyhow!("Failed to run nextflow: {}", e))?;

        if status.success() {
            println!("✅ Test completed successfully!");
            println!("   Results saved to: {}", output_dir.display());
        } else {
            println!("❌ Test failed with exit code: {:?}", status.code());
        }
    } else {
        println!("❌ No participants specified in submission");
    }

    Ok(())
}

// NOTE: Keep tests at the very end of file to satisfy clippy (items-after-test-module)
#[cfg(test)]
mod tests_final {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    fn write_submission(dir: &Path, sender: &str, filename: &str, status: &str) {
        let sender_dir = dir.join(sender);
        fs::create_dir_all(&sender_dir).unwrap();
        let path = sender_dir.join(filename);
        let yaml = format!(
            r#"name: Proj
author: A
datasites: ["d@example"]
participants: ["P1"]
syft_url: syft://x/y
status: {}
"#,
            status
        );
        fs::write(path, yaml).unwrap();
    }

    #[test]
    fn extract_date_from_filename_parses_and_falls_back() {
        assert_eq!(
            extract_date_from_filename("proj-2024-09-18-abc"),
            "2024-09-18"
        );
        assert_eq!(extract_date_from_filename("bad"), "unknown");
        assert_eq!(extract_date_from_filename("x-y-z"), "unknown");
    }

    #[test]
    fn get_inbox_path_uses_biovault_home() {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));
        let p = get_inbox_path().unwrap();
        assert!(p.ends_with("inbox"));
        assert!(p.starts_with(tmp.path().join(".bv")));
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn format_status_colors_known_statuses() {
        for s in ["pending", "approved", "rejected", "reviewing", "other"] {
            let out = format_status(s);
            assert!(!out.is_empty());
        }
    }

    #[test]
    fn load_submissions_filters_and_sorts() {
        let tmp = TempDir::new().unwrap();
        let home = tmp.path().join(".bv");
        crate::config::set_test_biovault_home(&home);
        let inbox = home.join("inbox");
        write_submission(
            &inbox,
            "alice@example.com",
            "proj-2024-09-18-aaaa.yaml",
            "pending",
        );
        write_submission(
            &inbox,
            "bob@example.com",
            "proj-2023-01-02-bbbb.yaml",
            "rejected",
        );
        let subs = load_submissions(false).unwrap();
        assert_eq!(subs.len(), 1);
        assert_eq!(subs[0].0, "alice@example.com");
        let subs_all = load_submissions(true).unwrap();
        assert_eq!(subs_all.len(), 2);
        assert_eq!(subs_all[0].0, "alice@example.com");
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn display_helpers_run_without_panic() {
        let sub = InboxSubmission {
            name: "Proj".into(),
            author: "A".into(),
            datasites: Some(vec!["d@example".into()]),
            participants: Some(vec!["P".into()]),
            syft_url: "syft://x/y".into(),
            status: "pending".into(),
        };
        let items = vec![(
            "sender@example".into(),
            sub.clone(),
            PathBuf::from("/tmp/id.yaml"),
            "2024-09-18".into(),
        )];
        display_concise_list(&items);
        display_full_submission(1, "sender@example", &sub, Path::new("/tmp/id.yaml"));
    }

    #[test]
    fn show_submission_detail_prints() {
        let sub = InboxSubmission {
            name: "Proj".into(),
            author: "A".into(),
            datasites: Some(vec!["d@example".into()]),
            participants: Some(vec![
                "P1".into(),
                "P2".into(),
                "P3".into(),
                "P4".into(),
                "P5".into(),
                "P6".into(),
            ]),
            syft_url: "syft://x/y".into(),
            status: "approved".into(),
        };
        let data = (
            "sender@example".into(),
            sub,
            PathBuf::from("/tmp/id.yaml"),
            "2024-09-18".into(),
        );
        show_submission_detail(&data);
    }

    #[tokio::test]
    async fn show_by_index_displays_detail() {
        let tmp = TempDir::new().unwrap();
        let home = tmp.path().join(".bv");
        crate::config::set_test_biovault_home(&home);
        let inbox = home.join("inbox");
        write_submission(
            &inbox,
            "alice@example.com",
            "proj-2024-09-18-aaaa.yaml",
            "pending",
        );
        super::show("1", true).await.unwrap();
        crate::config::clear_test_biovault_home();
    }

    #[test]
    fn show_submission_files_lists_without_error() {
        let tmp = TempDir::new().unwrap();
        // Build a fake submission directory structure
        let sub_dir = tmp
            .path()
            .join("sender@example.com")
            .join("proj-2024-01-01-deadbeef");
        std::fs::create_dir_all(sub_dir.join("assets/nested")).unwrap();
        std::fs::write(sub_dir.join("workflow.nf"), b"wf").unwrap();
        std::fs::write(sub_dir.join("assets/nested/file.txt"), b"x").unwrap();
        // The path stored in tuple points to a yaml file inside the submission dir
        let path = sub_dir.join("project.yaml");
        std::fs::write(&path, b"name: P\nauthor: A\nstatus: pending\nsyft_url: x\n").unwrap();
        let sub = InboxSubmission {
            name: "P".into(),
            author: "A".into(),
            datasites: None,
            participants: None,
            syft_url: "x".into(),
            status: "pending".into(),
        };
        let data = ("sender@example.com".into(), sub, path, "2024-01-01".into());
        show_submission_files(&data).unwrap();
    }

    #[tokio::test]
    async fn list_shows_concise_overview() {
        let tmp = TempDir::new().unwrap();
        let home = tmp.path().join(".bv");
        crate::config::set_test_biovault_home(&home);
        let inbox = home.join("inbox");
        write_submission(
            &inbox,
            "alice@example.com",
            "proj-2024-09-18-aaaa.yaml",
            "pending",
        );
        write_submission(
            &inbox,
            "bob@example.com",
            "proj-2023-01-02-bbbb.yaml",
            "approved",
        );
        super::list(false, false).await.unwrap();
        crate::config::clear_test_biovault_home();
    }

    #[tokio::test]
    async fn show_handles_invalid_index_and_name_matches() {
        let tmp = TempDir::new().unwrap();
        let home = tmp.path().join(".bv");
        crate::config::set_test_biovault_home(&home);
        let inbox = home.join("inbox");
        // Two submissions with same project name to trigger multi-match
        write_submission(
            &inbox,
            "alice@example.com",
            "proj-2024-09-18-aaaa.yaml",
            "pending",
        );
        write_submission(
            &inbox,
            "bob@example.com",
            "proj-2024-09-19-bbbb.yaml",
            "approved",
        );

        // Invalid index prints error but returns Ok
        super::show("5", true).await.unwrap();
        // Name-based no match
        super::show("does-not-exist", true).await.unwrap();
        // Name-based multi-match ('proj' present in both filenames)
        super::show("proj", true).await.unwrap();

        crate::config::clear_test_biovault_home();
    }

    #[tokio::test]
    async fn create_scaffold_project_without_example() {
        let tmp = TempDir::new().unwrap();
        let proj_dir = tmp.path().join("myproj");
        // Provide email via env var to avoid requiring a real config
        std::env::set_var("SYFTBOX_EMAIL", "scaffold@example.com");
        // Force non-interactive mode
        std::env::set_var("BIOVAULT_NON_INTERACTIVE", "1");

        super::create(
            Some("myproj".into()),
            Some(proj_dir.to_string_lossy().to_string()),
            None,
            None,
            None,
            None,
        )
        .await
        .unwrap();

        assert!(proj_dir.join("project.yaml").exists());
        assert!(proj_dir.join("workflow.nf").exists());
        assert!(proj_dir.join("assets").is_dir());

        std::env::remove_var("SYFTBOX_EMAIL");
        std::env::remove_var("BIOVAULT_NON_INTERACTIVE");
    }

    #[test]
    fn list_examples_runs() {
        // Should list embedded examples without error
        super::list_examples().unwrap();
    }

    #[test]
    fn test_format_status_pending() {
        let result = super::format_status("pending");
        assert!(result.contains("pending"));
    }

    #[test]
    fn test_format_status_approved() {
        let result = super::format_status("approved");
        assert!(result.contains("approved"));
    }

    #[test]
    fn test_format_status_rejected() {
        let result = super::format_status("rejected");
        assert!(result.contains("rejected"));
    }

    #[test]
    fn test_format_status_reviewing() {
        let result = super::format_status("reviewing");
        assert!(result.contains("reviewing"));
    }

    #[test]
    fn test_format_status_unknown() {
        let result = super::format_status("unknown_status");
        assert_eq!(result, "unknown_status");
    }

    #[test]
    fn test_get_inbox_path() {
        let tmp = TempDir::new().unwrap();
        crate::config::set_test_biovault_home(tmp.path().join(".bv"));

        let result = super::get_inbox_path();
        assert!(result.is_ok());
        let path = result.unwrap();
        // Path should contain inbox
        assert!(path.to_string_lossy().contains("inbox"));

        crate::config::clear_test_biovault_home();
    }

    #[tokio::test]
    async fn test_show_empty_submissions() {
        let result = super::show("1", false).await;
        // May error or succeed depending on state - just verify it doesn't panic
        let _ = result;
    }
}
