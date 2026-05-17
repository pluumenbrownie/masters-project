use mcm_finder_lib::dataset::VecDataset;

use crate::app::App;

pub mod app;
pub mod event;
pub mod ui;

#[tokio::main]
async fn main() -> color_eyre::Result<()> {
    color_eyre::install()?;
    let terminal = ratatui::init();
    let result = App::<VecDataset>::new().run(terminal).await;
    ratatui::restore();
    result
}
