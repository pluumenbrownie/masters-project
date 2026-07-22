use std::sync::mpsc::{self, Sender};

use crate::collector::Collector;
use crate::event::CollectorEvent;
use crate::handler::{FnHandler, Handler};

/// Configures and builds a [`Collector`] together with its [`Sender`].
///
/// # Example
/// ```rust
/// use annolog::{CollectorBuilder, CollectorEvent};
///
/// let (collector, tx) = CollectorBuilder::<String>::new()
///     .with_fn(|s| println!("{s}"))
///     .build();
///
/// std::thread::spawn(|| collector.run());
/// tx.send(CollectorEvent::Data("hello".into())).unwrap();
/// tx.send(CollectorEvent::Shutdown).unwrap();
/// ```
pub struct CollectorBuilder<T> {
    handlers: Vec<Box<dyn Handler<T>>>,
}

impl<T: Send + 'static> CollectorBuilder<T> {
    /// Creates an empty builder with no handlers.
    pub fn new() -> Self {
        Self {
            handlers: Vec::new(),
        }
    }

    /// Registers a handler that implements [`Handler<T>`].
    pub fn with_handler(mut self, handler: impl Handler<T> + 'static) -> Self {
        self.handlers.push(Box::new(handler));
        self
    }

    /// Registers a closure as a handler — no [`Handler`] impl required.
    pub fn with_fn(self, f: impl FnMut(&T) + Send + 'static) -> Self {
        self.with_handler(FnHandler::new(f))
    }

    /// Consumes the builder and returns a `(Collector<T>, Sender<CollectorEvent<T>>)` pair.
    ///
    /// - Keep the [`Sender`] (and clone it freely) to send events from producers.
    /// - Pass the [`Collector`] to a dedicated thread and call `.run()`.
    pub fn build(self) -> (Collector<T>, Sender<CollectorEvent<T>>) {
        let (tx, rx) = mpsc::channel();
        let collector = Collector {
            receiver: rx,
            handlers: self.handlers,
        };
        (collector, tx)
    }
}

impl<T: Send + 'static> Default for CollectorBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}
