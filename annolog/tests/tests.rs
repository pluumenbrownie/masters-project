use std::sync::{Arc, Mutex};
use std::time::Duration;

use annolog::{BufferingHandler, CollectorBuilder, CollectorEvent, FnHandler, Handler};

// ---------------------------------------------------------------------------
// User-defined types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
struct UserStruct {
    id: u32,
    label: String,
}

#[derive(Debug, Clone, PartialEq)]
enum MetricEvent {
    Counter { name: String, value: u64 },
    Gauge { name: String, value: f64 },
}

#[derive(Debug, Clone, PartialEq)]
enum LogEvent {
    Info(String),
    Error { message: String, code: u32 },
}

#[derive(Debug, Clone, PartialEq)]
enum AppEvent {
    Metric(MetricEvent),
    Log(LogEvent),
    Raw(UserStruct),
}

// ---------------------------------------------------------------------------
// Recording handler
// ---------------------------------------------------------------------------

struct RecordingHandler {
    received: Arc<Mutex<Vec<AppEvent>>>,
}

impl RecordingHandler {
    fn new(store: Arc<Mutex<Vec<AppEvent>>>) -> Self {
        Self { received: store }
    }
}

impl Handler<AppEvent> for RecordingHandler {
    fn handle(&mut self, event: &AppEvent) {
        self.received.lock().unwrap().push(event.clone());
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn wait_for(store: &Arc<Mutex<Vec<AppEvent>>>, count: usize) {
    wait_for_timeout(store, count, Duration::from_secs(2));
}

fn wait_for_timeout(store: &Arc<Mutex<Vec<AppEvent>>>, count: usize, timeout: Duration) {
    let deadline = std::time::Instant::now() + timeout;
    loop {
        if store.lock().unwrap().len() >= count {
            break;
        }
        assert!(
            std::time::Instant::now() < deadline,
            "timed out waiting for {count} events"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
}

// ---------------------------------------------------------------------------
// Original tests
// ---------------------------------------------------------------------------

#[test]
fn test_single_event_received() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "hello".into(),
    ))))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store, 1);
    let received = store.lock().unwrap();
    assert_eq!(received.len(), 1);
    assert_eq!(received[0], AppEvent::Log(LogEvent::Info("hello".into())));
}

#[test]
fn test_multiple_producers() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    let tx2 = tx.clone();
    let tx3 = tx.clone();

    let h1 = std::thread::spawn(move || {
        tx2.send(CollectorEvent::Data(AppEvent::Metric(
            MetricEvent::Counter {
                name: "requests".into(),
                value: 1,
            },
        )))
        .unwrap();
    });

    let h2 = std::thread::spawn(move || {
        tx3.send(CollectorEvent::Data(AppEvent::Metric(MetricEvent::Gauge {
            name: "cpu".into(),
            value: 0.75,
        })))
        .unwrap();
    });

    tx.send(CollectorEvent::Data(AppEvent::Raw(UserStruct {
        id: 42,
        label: "main".into(),
    })))
    .unwrap();

    h1.join().unwrap();
    h2.join().unwrap();

    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store, 3);
    let received = store.lock().unwrap();
    assert_eq!(received.len(), 3);
    assert!(
        received
            .iter()
            .any(|e| matches!(e, AppEvent::Metric(MetricEvent::Counter { .. })))
    );
    assert!(
        received
            .iter()
            .any(|e| matches!(e, AppEvent::Metric(MetricEvent::Gauge { .. })))
    );
    assert!(received.iter().any(|e| matches!(e, AppEvent::Raw(_))));
}

#[test]
fn test_nested_enum_dispatch() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Error {
        message: "disk full".into(),
        code: 500,
    })))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store, 1);
    let received = store.lock().unwrap();

    match &received[0] {
        AppEvent::Log(LogEvent::Error { message, code }) => {
            assert_eq!(message, "disk full");
            assert_eq!(*code, 500);
        }
        other => panic!("unexpected event: {other:?}"),
    }
}

#[test]
fn test_closure_handler() {
    let store: Arc<Mutex<Vec<AppEvent>>> = Arc::new(Mutex::new(vec![]));
    let store_clone = store.clone();

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_fn(move |event| {
            store_clone.lock().unwrap().push(event.clone());
        })
        .build();

    std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "from closure".into(),
    ))))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store, 1);
    let received = store.lock().unwrap();
    assert_eq!(
        received[0],
        AppEvent::Log(LogEvent::Info("from closure".into()))
    );
}

#[test]
fn test_fanout_to_multiple_handlers() {
    let store1 = Arc::new(Mutex::new(vec![]));
    let store2 = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store1.clone()))
        .with_handler(RecordingHandler::new(store2.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Metric(
        MetricEvent::Counter {
            name: "hits".into(),
            value: 99,
        },
    )))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store1, 1);
    wait_for(&store2, 1);

    assert_eq!(store1.lock().unwrap().len(), 1);
    assert_eq!(store2.lock().unwrap().len(), 1);
    assert_eq!(store1.lock().unwrap()[0], store2.lock().unwrap()[0]);
}

#[test]
fn test_shutdown_stops_collector() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    let handle = std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "before shutdown".into(),
    ))))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();

    handle.join().expect("collector thread panicked");

    let _ = tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "after shutdown".into(),
    ))));

    let received = store.lock().unwrap();
    assert_eq!(received.len(), 1);
    assert_eq!(
        received[0],
        AppEvent::Log(LogEvent::Info("before shutdown".into()))
    );
}

// ---------------------------------------------------------------------------
// New tests
// ---------------------------------------------------------------------------

/// Dropping all senders without sending Shutdown should also stop the collector.
#[test]
fn test_channel_close_without_shutdown() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    let handle = std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "only event".into(),
    ))))
    .unwrap();

    // Drop the only sender — channel closes, collector should exit
    drop(tx);

    handle.join().expect("collector thread panicked");

    let received = store.lock().unwrap();
    assert_eq!(received.len(), 1);
    assert_eq!(
        received[0],
        AppEvent::Log(LogEvent::Info("only event".into()))
    );
}

/// Data events queued after Shutdown in the channel must not be processed.
#[test]
fn test_data_after_shutdown_not_processed() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    // Queue everything before spawning so ordering is deterministic
    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "before".into(),
    ))))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();
    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "after".into(),
    ))))
    .unwrap();

    let handle = std::thread::spawn(|| collector.run());
    handle.join().expect("collector thread panicked");

    let received = store.lock().unwrap();
    assert_eq!(received.len(), 1);
    assert_eq!(received[0], AppEvent::Log(LogEvent::Info("before".into())));
}

/// A handler that mutates its own state across events.
#[test]
fn test_stateful_handler() {
    struct CountingHandler {
        count: Arc<Mutex<u64>>,
    }

    impl Handler<AppEvent> for CountingHandler {
        fn handle(&mut self, event: &AppEvent) {
            if let AppEvent::Metric(MetricEvent::Counter { value, .. }) = event {
                *self.count.lock().unwrap() += value;
            }
        }
    }

    let total = Arc::new(Mutex::new(0u64));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(CountingHandler {
            count: total.clone(),
        })
        .build();

    let handle = std::thread::spawn(|| collector.run());

    for i in 1..=5 {
        tx.send(CollectorEvent::Data(AppEvent::Metric(
            MetricEvent::Counter {
                name: "x".into(),
                value: i,
            },
        )))
        .unwrap();
    }
    tx.send(CollectorEvent::Shutdown).unwrap();
    handle.join().unwrap();

    // 1+2+3+4+5 = 15
    assert_eq!(*total.lock().unwrap(), 15);
}

/// BufferingHandler with fewer events than capacity — partial batch is flushed on drop.
/// NOTE: this test documents current behaviour. If it fails, BufferingHandler
/// needs a Drop impl or the collector needs to flush handlers after the loop.
#[test]
fn test_buffering_handler_partial_flush() {
    let batches: Arc<Mutex<Vec<Vec<AppEvent>>>> = Arc::new(Mutex::new(vec![]));
    let batches_clone = batches.clone();

    let buffer = BufferingHandler::new(
        10, // capacity larger than number of events we'll send
        FnHandler::new(move |batch: &Vec<AppEvent>| {
            batches_clone.lock().unwrap().push(batch.clone());
        }),
    );

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(buffer)
        .build();

    let handle = std::thread::spawn(|| collector.run());

    // Send only 3 events, never reaching the capacity of 10
    for i in 0..3 {
        tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
            format!("msg-{i}"),
        ))))
        .unwrap();
    }
    tx.send(CollectorEvent::Shutdown).unwrap();
    handle.join().unwrap();

    let received = batches.lock().unwrap();
    assert_eq!(
        received.len(),
        1,
        "partial buffer should be flushed on shutdown"
    );
    assert_eq!(received[0].len(), 3);
}

/// BufferingHandler flushes exactly when capacity is reached, then again for remainder.
#[test]
fn test_buffering_handler_multiple_flushes() {
    let batches: Arc<Mutex<Vec<Vec<AppEvent>>>> = Arc::new(Mutex::new(vec![]));
    let batches_clone = batches.clone();

    let buffer = BufferingHandler::new(
        3,
        FnHandler::new(move |batch: &Vec<AppEvent>| {
            batches_clone.lock().unwrap().push(batch.clone());
        }),
    );

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(buffer)
        .build();

    let handle = std::thread::spawn(|| collector.run());

    // 7 events: flush at 3, flush at 6, then 1 remainder flushed on shutdown
    for i in 0..7 {
        tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
            format!("msg-{i}"),
        ))))
        .unwrap();
    }
    tx.send(CollectorEvent::Shutdown).unwrap();
    handle.join().unwrap();

    let received = batches.lock().unwrap();
    assert_eq!(received.len(), 3);
    assert_eq!(received[0].len(), 3);
    assert_eq!(received[1].len(), 3);
    assert_eq!(received[2].len(), 1);
}

/// High-volume: 10,000 events from 10 threads all arrive with no drops.
#[test]
fn test_stress_high_volume() {
    const THREADS: usize = 10;
    const PER_THREAD: usize = 1_000;
    const TOTAL: usize = THREADS * PER_THREAD;

    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    let handles: Vec<_> = (0..THREADS)
        .map(|t| {
            let tx = tx.clone();
            std::thread::spawn(move || {
                for i in 0..PER_THREAD {
                    tx.send(CollectorEvent::Data(AppEvent::Metric(
                        MetricEvent::Counter {
                            name: format!("thread-{t}"),
                            value: i as u64,
                        },
                    )))
                    .unwrap();
                }
            })
        })
        .collect();

    for h in handles {
        h.join().unwrap();
    }
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for_timeout(&store, TOTAL, Duration::from_secs(10));
    assert_eq!(store.lock().unwrap().len(), TOTAL);
}

/// 50 threads each send exactly 1 event — all 50 must arrive.
#[test]
fn test_many_senders_one_event_each() {
    const SENDERS: usize = 50;

    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    let handles: Vec<_> = (0..SENDERS)
        .map(|i| {
            let tx = tx.clone();
            std::thread::spawn(move || {
                tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
                    format!("sender-{i}"),
                ))))
                .unwrap();
            })
        })
        .collect();

    for h in handles {
        h.join().unwrap();
    }
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store, SENDERS);
    assert_eq!(store.lock().unwrap().len(), SENDERS);
}

/// Building with no handlers should not panic — events are silently consumed.
#[test]
fn test_no_handlers_does_not_panic() {
    let (collector, tx) = CollectorBuilder::<AppEvent>::new().build();

    let handle = std::thread::spawn(|| collector.run());

    for i in 0..10 {
        tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
            format!("event-{i}"),
        ))))
        .unwrap();
    }
    tx.send(CollectorEvent::Shutdown).unwrap();

    handle.join().expect("collector panicked with no handlers");
}

/// CollectorBuilder::default() is equivalent to CollectorBuilder::new().
#[test]
fn test_builder_default_equivalent_to_new() {
    let store = Arc::new(Mutex::new(vec![]));

    let (collector, tx) = CollectorBuilder::<AppEvent>::default()
        .with_handler(RecordingHandler::new(store.clone()))
        .build();

    std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "via default".into(),
    ))))
    .unwrap();
    tx.send(CollectorEvent::Shutdown).unwrap();

    wait_for(&store, 1);
    assert_eq!(
        store.lock().unwrap()[0],
        AppEvent::Log(LogEvent::Info("via default".into()))
    );
}

/// A panicking handler causes the collector thread to panic, which can be
/// observed via JoinHandle::join returning Err.
#[test]
fn test_panicking_handler_propagates() {
    struct PanicHandler;

    impl Handler<AppEvent> for PanicHandler {
        fn handle(&mut self, _event: &AppEvent) {
            panic!("handler intentionally panicked");
        }
    }

    let (collector, tx) = CollectorBuilder::<AppEvent>::new()
        .with_handler(PanicHandler)
        .build();

    let handle = std::thread::spawn(|| collector.run());

    tx.send(CollectorEvent::Data(AppEvent::Log(LogEvent::Info(
        "trigger".into(),
    ))))
    .unwrap();

    // The collector thread should have panicked
    let result = handle.join();
    assert!(result.is_err(), "expected collector thread to panic");
}
