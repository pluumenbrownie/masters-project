use std::sync::mpsc::Receiver;

use crate::event::CollectorEvent;
use crate::handler::Handler;

/// The collector owns the receiving end of the MPSC channel and a list of
/// handlers. Call [`Collector::run`] on a dedicated thread to start processing.
pub struct Collector<T> {
    pub(crate) receiver: Receiver<CollectorEvent<T>>,
    pub(crate) handlers: Vec<Box<dyn Handler<T>>>,
}

impl<T: Send + 'static> Collector<T> {
    /// Blocking event loop. Exits when a [`CollectorEvent::Shutdown`] is
    /// received or all senders are dropped (channel closed).
    pub fn run(mut self) {
        for event in self.receiver.iter() {
            match event {
                CollectorEvent::Shutdown => break,
                CollectorEvent::Data(data) => {
                    for handler in &mut self.handlers {
                        handler.handle(&data);
                    }
                }
            }
        }
    }
}
