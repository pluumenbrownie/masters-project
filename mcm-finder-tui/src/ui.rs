use mcm_finder_lib::dataset::Dataset;
use ratatui::{
    buffer::Buffer,
    layout::{Alignment, Constraint, Layout, Rect},
    style::{Color, Stylize},
    widgets::{Block, BorderType, List, Paragraph, StatefulWidget, Widget},
};

use crate::app::{App, AppState};

impl<T: Dataset> Widget for &mut App<T> {
    /// Renders the user interface widgets.
    ///
    fn render(self, area: Rect, buf: &mut Buffer) {
        // let layout = Layout::default().direction(Direction::Vertical).split(area);
        let block = Block::bordered()
            .title("<MCM finder TUI>")
            .title_alignment(Alignment::Center)
            .border_type(BorderType::Rounded);
        let block_area = block.inner(area);
        block.render(area, buf);

        self.render_inner(block_area, buf);
    }
}

impl<T: Dataset> App<T> {
    fn render_inner(&mut self, area: Rect, buf: &mut Buffer) {
        match self.state {
            AppState::Default { mut list_state } => {
                let inner_block = Block::bordered()
                    .title("<Dataset>")
                    .title_alignment(Alignment::Center)
                    .border_type(BorderType::Rounded);
                inner_block.render(area, buf);

                if let Some(_dataset) = &self.dataset {
                    todo!()
                } else {
                    let text = format!("No dataset loaded. {:?}", list_state);
                    let paragraph = Paragraph::new(text).fg(Color::Cyan);

                    let list = List::new(self.get_file_list()).highlight_symbol(">");

                    let layout = Layout::default()
                        .constraints([Constraint::Length(1), Constraint::Length(list.len() as u16)])
                        .spacing(1)
                        .margin(1)
                        .split(area);

                    paragraph.render(layout[0], buf);
                    StatefulWidget::render(list, layout[1], buf, &mut list_state);
                };
            }
        }
    }

    fn get_file_list(&self) -> Vec<String> {
        vec![
            String::from("./mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat"),
            String::from("./mcm-finder-lib/tests/data/MNIST11.sorted"),
            String::from("./mcm-finder-lib/tests/data/Big5PT.sorted"),
        ]
    }
}
