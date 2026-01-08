use opentelemetry::global;
use opentelemetry::trace::TracerProvider as _;
use opentelemetry::KeyValue;
use opentelemetry_otlp::WithExportConfig;
use opentelemetry_sdk::trace::TracerProvider;
use std::env;
use std::sync::OnceLock;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

static TRACER_PROVIDER: OnceLock<Option<TracerProvider>> = OnceLock::new();
static TELEMETRY_RUNTIME: OnceLock<tokio::runtime::Runtime> = OnceLock::new();

pub fn init() {
    TRACER_PROVIDER.get_or_init(|| {
        // Telemetry is opt-in: requires both OTEL_EXPORTER_OTLP_ENDPOINT and BIOVAULT_ENABLE_TELEMETRY=1
        if env::var("BIOVAULT_ENABLE_TELEMETRY").unwrap_or_default() != "1" {
            return None;
        }

        let endpoint = match env::var("OTEL_EXPORTER_OTLP_ENDPOINT") {
            Ok(e) if !e.is_empty() => e,
            _ => return None,
        };

        let service_name =
            env::var("OTEL_SERVICE_NAME").unwrap_or_else(|_| "biovault-cli".to_string());

        match try_init_tracer(&endpoint, &service_name) {
            Ok(provider) => {
                tracing::info!(
                    "Telemetry enabled: endpoint={}, service={}",
                    endpoint,
                    service_name
                );
                Some(provider)
            }
            Err(e) => {
                tracing::warn!("Failed to initialize telemetry: {}", e);
                None
            }
        }
    });
}

fn try_init_tracer(
    endpoint: &str,
    service_name: &str,
) -> Result<TracerProvider, Box<dyn std::error::Error>> {
    // Create a dedicated runtime for telemetry batch export
    let rt = TELEMETRY_RUNTIME.get_or_init(|| {
        tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .thread_name("otel-export")
            .enable_all()
            .build()
            .expect("Failed to create telemetry runtime")
    });

    // Enter the runtime context to allow batch exporter creation
    let _guard = rt.enter();

    let exporter = opentelemetry_otlp::SpanExporter::builder()
        .with_http()
        .with_endpoint(format!("{}/v1/traces", endpoint))
        .build()?;

    let provider = TracerProvider::builder()
        .with_batch_exporter(exporter, opentelemetry_sdk::runtime::Tokio)
        .with_resource(opentelemetry_sdk::Resource::new(vec![KeyValue::new(
            "service.name",
            service_name.to_string(),
        )]))
        .build();

    global::set_tracer_provider(provider.clone());

    let telemetry_layer =
        tracing_opentelemetry::layer().with_tracer(provider.tracer("biovault-cli"));

    tracing_subscriber::registry()
        .with(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .with(telemetry_layer)
        .try_init()
        .ok();

    Ok(provider)
}

#[allow(dead_code)]
pub fn shutdown() {
    if let Some(Some(provider)) = TRACER_PROVIDER.get() {
        if let Err(e) = provider.shutdown() {
            eprintln!("Error shutting down tracer: {:?}", e);
        }
    }
}
