use std::ffi::OsString;

use crate::event::{AppEvent, Event, EventHandler};
use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use mcm_finder_lib::dataset::Dataset;
use ratatui::{DefaultTerminal, widgets::ListState};

/// Application.
#[derive(Debug)]
pub struct App<T: Dataset> {
    /// Is the application running?
    pub running: bool,
    /// Event handler.
    pub events: EventHandler,
    /// Current app state.
    pub state: AppState,
    /// The path to the dataset
    pub data_path: Option<DatasetPath>,

    pub dataset: Option<T>,
}

impl<T: Dataset> Default for App<T> {
    fn default() -> Self {
        Self {
            running: true,
            events: EventHandler::new(),
            state: AppState::default(),
            data_path: None,
            dataset: None,
        }
    }
}

impl<T: Dataset> App<T> {
    /// Constructs a new instance of [`App`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Run the application's main loop.
    pub async fn run(mut self, mut terminal: DefaultTerminal) -> color_eyre::Result<()> {
        while self.running {
            terminal.draw(|frame| frame.render_widget(&mut self, frame.area()))?;
            match self.events.next().await? {
                Event::Tick => self.tick(),
                Event::Crossterm(event) => match event {
                    crossterm::event::Event::Key(key_event)
                        if key_event.kind == crossterm::event::KeyEventKind::Press =>
                    {
                        self.handle_key_events(key_event)?
                    }
                    _ => {}
                },
                Event::App(app_event) => match app_event {
                    AppEvent::Quit => self.quit(),
                    AppEvent::ScrollUp => self.up(),
                    AppEvent::ScrollDown => self.down(),
                },
            }
        }
        Ok(())
    }

    /// Handles the key events and updates the state of [`App`].
    pub fn handle_key_events(&mut self, key_event: KeyEvent) -> color_eyre::Result<()> {
        match key_event.code {
            KeyCode::Esc | KeyCode::Char('q') => self.events.send(AppEvent::Quit),
            KeyCode::Char('c' | 'C') if key_event.modifiers == KeyModifiers::CONTROL => {
                self.events.send(AppEvent::Quit)
            }

            KeyCode::Up | KeyCode::Char('k') => self.events.send(AppEvent::ScrollUp),
            KeyCode::Down | KeyCode::Char('j') => self.events.send(AppEvent::ScrollDown),
            // Other handlers you could add here.
            _ => {}
        }
        Ok(())
    }

    /// Handles the tick event of the terminal.
    ///
    /// The tick event is where you can update the state of your application with any logic that
    /// needs to be updated at a fixed frame rate. E.g. polling a server, updating an animation.
    pub fn tick(&self) {}

    /// Set running to false to quit the application.
    pub fn quit(&mut self) {
        self.running = false;
    }

    fn up(&mut self) {
        match self.state {
            AppState::Default { mut list_state } => list_state.select_previous(),
        }
    }

    fn down(&mut self) {
        match self.state {
            AppState::Default { mut list_state } => list_state.select_next(),
        }
    }
}

#[derive(Debug)]
pub enum AppState {
    Default { list_state: ListState },
}

impl Default for AppState {
    fn default() -> Self {
        AppState::Default {
            list_state: ListState::default().with_selected(Some(0)),
        }
    }
}

#[derive(Debug)]
pub enum DatasetPath {
    SupremeCourt,
    MnistTruncated,
    Big5,
    Custom(OsString),
}
