//! Recomputed full phase-domain operators for the CoupledGrid example.
use openline::{Catalog, Line, Phase, Transposition};
use openline_compute::{
    CsvWriter, FormulationSet, FrequencySweep, LineStudy, MatrixDomain,
    MatrixQuantityKind, Model, OperatingPoint,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output = std::env::args().nth(1).expect("output CSV path required");
    let catalog = Catalog::bundled()?;
    let line = Line::builder("coupled_13_8kv")
        .description("Untransposed three-wire ACSR Linnet distribution feeder")
        .cross_section("linnet", |xs| {
            xs.circuit("circuit_1", |c| {
                c.phase_m(Phase::A, "phase_a", -1.2, 10.5)
                    .phase_m(Phase::B, "phase_b", 0.0, 11.0)
                    .phase_m(Phase::C, "phase_c", 1.2, 10.5)
                    .conductor("ACSR_Linnet")
            })
            .transposition(Transposition::Untransposed)
        })
        .segment_m("reference", "linnet", 1000.0)
        .build(&catalog)?;
    let result = LineStudy::new(&line)
        .segment("reference")
        .operating_point(
            OperatingPoint::overhead()
                .conductor_temperature_c(20.0)
                .ambient_temperature_c(20.0)
                .earth_resistivity_ohm_m(100.0)
                .with_electrical(|e| e.system_voltage_kv(13.8)),
        )
        .frequencies(FrequencySweep::log_hz(1.0, 30000.0, 801)?)
        .formulations(FormulationSet::new(vec![
            Model::SchelkunoffInternal,
            Model::PerfectEarthExternal,
            Model::CarsonEarthReturn,
            Model::MaxwellPerfectEarth,
        ]))
        .matrix_domain(MatrixDomain::Phase)
        .run()?;
    CsvWriter::new(vec![
        MatrixQuantityKind::SeriesImpedance,
        MatrixQuantityKind::ShuntAdmittance,
    ])
    .write_path(&result, output)?;
    Ok(())
}
