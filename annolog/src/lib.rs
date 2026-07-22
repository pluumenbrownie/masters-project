//! # annolog
//!
//! A generic, handler-based event collector over MPSC channels.
//!
//! `annolog` lets you send any of your own Rust types — structs, enums, nested
//! enums — from anywhere in your codebase to a central [`Collector`] running on
//! a dedicated thread. One or more [`Handler`]s registered on the collector
//! decide what to do with each event.
//!
//! ## Core concepts
//!
//! | Type | Role |
//! |---|---|
//! | [`CollectorEvent<T>`] | Thin envelope around your type `T`. Carries either `Data(T)` or a `Shutdown` signal. |
//! | [`Collector<T>`] | Owns the receiving end of the channel and the list of handlers. Call `.run()` on a dedicated thread. |
//! | [`CollectorBuilder<T>`] | Fluent builder that wires up handlers and produces a `(Collector<T>, Sender<CollectorEvent<T>>)` pair. |
//! | [`Handler<T>`] | Trait you implement (or use a built-in) to process incoming events. |
//!
//! ## Quick start
//!
//! ```rust
//! use annolog::{CollectorBuilder, CollectorEvent};
//!
//! let (collector, tx) = CollectorBuilder::<String>::new()
//!     .with_fn(|s| println!("{s}"))
//!     .build();
//!
//! std::thread::spawn(|| collector.run());
//!
//! tx.send(CollectorEvent::Data("hello".into())).unwrap();
//! tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! ## Defining your own event types
//!
//! `annolog` is generic over `T` — your type just needs to be `Send + 'static`.
//! Use plain structs, enums, or deeply nested enums:
//!
//! ```rust
//! #[derive(Debug, Clone)]
//! pub struct RequestInfo {
//!     pub path: String,
//!     pub status: u16,
//! }
//!
//! #[derive(Debug, Clone)]
//! pub enum MetricEvent {
//!     Counter { name: String, value: u64 },
//!     Gauge   { name: String, value: f64 },
//! }
//!
//! #[derive(Debug, Clone)]
//! pub enum AppEvent {
//!     Metric(MetricEvent),   // nested enum
//!     Request(RequestInfo),  // plain struct
//!     Shutdown,              // your own domain signal (distinct from CollectorEvent::Shutdown)
//! }
//! ```
//!
//! ## Implementing a custom handler
//!
//! Implement [`Handler<T>`] for any type that should process events.
//! The handler runs on the collector thread and has exclusive `&mut self` access,
//! so it can freely hold and mutate local state without any locking:
//!
//! ```rust
//! use annolog::Handler;
//! # #[derive(Debug, Clone)] enum AppEvent { Count(u64) }
//!
//! struct SummingHandler {
//!     total: u64,
//! }
//!
//! impl Handler<AppEvent> for SummingHandler {
//!     fn handle(&mut self, event: &AppEvent) {
//!         if let AppEvent::Count(n) = event {
//!             self.total += n;
//!         }
//!     }
//! }
//! ```
//!
//! ## Using a closure instead
//!
//! For simple cases, skip the `Handler` impl entirely and pass a closure via
//! [`CollectorBuilder::with_fn`]:
//!
//! ```rust
//! use annolog::{CollectorBuilder, CollectorEvent};
//!
//! let (collector, tx) = CollectorBuilder::<String>::new()
//!     .with_fn(|event| eprintln!("[log] {event}"))
//!     .build();
//! # std::thread::spawn(|| collector.run());
//! # tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! ## Multiple handlers (fanout)
//!
//! Chain as many handlers as you need. Every handler receives every event in
//! registration order:
//!
//! ```rust
//! use annolog::{CollectorBuilder, CollectorEvent};
//! # use annolog::Handler;
//! # #[derive(Debug, Clone)] enum AppEvent { Foo }
//! # struct DbHandler; impl Handler<AppEvent> for DbHandler { fn handle(&mut self, _: &AppEvent) {} }
//! # struct MetricsHandler; impl Handler<AppEvent> for MetricsHandler { fn handle(&mut self, _: &AppEvent) {} }
//!
//! let (collector, tx) = CollectorBuilder::<AppEvent>::new()
//!     .with_handler(DbHandler)
//!     .with_handler(MetricsHandler)
//!     .with_fn(|e| println!("{e:?}"))
//!     .build();
//! # std::thread::spawn(|| collector.run());
//! # tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! ## Multiple producers
//!
//! [`Sender<CollectorEvent<T>>`] is `Clone` — hand out copies to as many
//! producer threads as you need:
//!
//! ```rust
//! use annolog::{CollectorBuilder, CollectorEvent};
//!
//! let (collector, tx) = CollectorBuilder::<String>::new()
//!     .with_fn(|s| println!("{s}"))
//!     .build();
//!
//! std::thread::spawn(|| collector.run());
//!
//! let handles: Vec<_> = (0..4).map(|i| {
//!     let tx = tx.clone();
//!     std::thread::spawn(move || {
//!         tx.send(CollectorEvent::Data(format!("from thread {i}"))).unwrap();
//!     })
//! }).collect();
//!
//! for h in handles { h.join().unwrap(); }
//! tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! ## Filtering events
//!
//! Wrap any handler in [`FilterHandler`] to only forward matching events:
//!
//! ```rust
//! use annolog::{CollectorBuilder, CollectorEvent, FilterHandler};
//! # use annolog::Handler;
//! # #[derive(Debug, Clone)] enum AppEvent { Log(String), Metric(u64) }
//! # struct LogHandler; impl Handler<AppEvent> for LogHandler { fn handle(&mut self, _: &AppEvent) {} }
//!
//! let filter = FilterHandler::new(
//!     |e| matches!(e, AppEvent::Log(_)),  // predicate
//!     LogHandler,                          // inner handler, only called when predicate is true
//! );
//!
//! let (collector, tx) = CollectorBuilder::<AppEvent>::new()
//!     .with_handler(filter)
//!     .build();
//! # std::thread::spawn(|| collector.run());
//! # tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! ## Batching events
//!
//! Wrap any `Handler<Vec<T>>` in [`BufferingHandler`] to accumulate events and
//! forward them in batches. Any events still in the buffer when the collector
//! shuts down are automatically flushed:
//!
//! ```rust
//! use annolog::{BufferingHandler, CollectorBuilder, CollectorEvent, FnHandler};
//!
//! let buffer = BufferingHandler::new(
//!     100,  // flush every 100 events
//!     FnHandler::new(|batch: &Vec<String>| {
//!         println!("flushing {} events", batch.len());
//!     }),
//! );
//!
//! let (collector, tx) = CollectorBuilder::<String>::new()
//!     .with_handler(buffer)
//!     .build();
//! # std::thread::spawn(|| collector.run());
//! # tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! ## Stopping the collector
//!
//! Send [`CollectorEvent::Shutdown`] to stop the collector after it has
//! processed all events queued before the signal:
//!
//! ```rust
//! # use annolog::{CollectorBuilder, CollectorEvent};
//! # let (collector, tx) = CollectorBuilder::<String>::new().build();
//! # std::thread::spawn(|| collector.run());
//! tx.send(CollectorEvent::Shutdown).unwrap();
//! ```
//!
//! Alternatively, drop all [`Sender`] handles — the collector exits cleanly
//! once the channel is closed.
//!
//! [`Sender`]: std::sync::mpsc::Sender

mod builder;
mod collector;
mod event;
mod handler;

pub use builder::CollectorBuilder;
pub use collector::Collector;
pub use event::CollectorEvent;
pub use handler::{BufferingHandler, FilterHandler, FnHandler, Handler};
