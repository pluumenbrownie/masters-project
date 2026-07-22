/// Trait implemented by anything that can receive and process a collected event.
///
/// # Example
/// ```rust
/// use annolog::Handler;
///
/// struct PrintHandler;
///
/// impl Handler<String> for PrintHandler {
///     fn handle(&mut self, event: &String) {
///         println!("{event}");
///     }
/// }
/// ```
pub trait Handler<T>: Send + 'static {
    fn handle(&mut self, event: &T);
}

// ---------------------------------------------------------------------------
// FnHandler — wraps a closure so users don't need a full Handler impl
// ---------------------------------------------------------------------------

pub struct FnHandler<T> {
    f: Box<dyn FnMut(&T) + Send + 'static>,
}

impl<T> FnHandler<T> {
    pub fn new(f: impl FnMut(&T) + Send + 'static) -> Self {
        Self { f: Box::new(f) }
    }
}

impl<T: 'static> Handler<T> for FnHandler<T> {
    fn handle(&mut self, event: &T) {
        (self.f)(event);
    }
}

// ---------------------------------------------------------------------------
// FilterHandler — only forwards events that pass a predicate
// ---------------------------------------------------------------------------

pub struct FilterHandler<T, H: Handler<T>> {
    predicate: Box<dyn Fn(&T) -> bool + Send + 'static>,
    inner: H,
    _marker: std::marker::PhantomData<T>,
}

impl<T, H: Handler<T>> FilterHandler<T, H> {
    pub fn new(predicate: impl Fn(&T) -> bool + Send + 'static, inner: H) -> Self {
        Self {
            predicate: Box::new(predicate),
            inner,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<T: Send + 'static, H: Handler<T>> Handler<T> for FilterHandler<T, H> {
    fn handle(&mut self, event: &T) {
        if (self.predicate)(event) {
            self.inner.handle(event);
        }
    }
}

// ---------------------------------------------------------------------------
// BufferingHandler — batches N events before forwarding the batch.
// Any remaining events in the buffer are flushed when the handler is dropped,
// which happens when the Collector exits its run loop.
// ---------------------------------------------------------------------------

pub struct BufferingHandler<T: Clone, H: Handler<Vec<T>>> {
    buffer: Vec<T>,
    capacity: usize,
    inner: H,
}

impl<T: Clone, H: Handler<Vec<T>>> BufferingHandler<T, H> {
    pub fn new(capacity: usize, inner: H) -> Self {
        assert!(capacity > 0, "BufferingHandler capacity must be > 0");
        Self {
            buffer: Vec::with_capacity(capacity),
            capacity,
            inner,
        }
    }

    /// Flush whatever is in the buffer, even if not full.
    pub fn flush(&mut self) {
        if !self.buffer.is_empty() {
            let batch = std::mem::take(&mut self.buffer);
            self.inner.handle(&batch);
        }
    }
}

impl<T: Clone + Send + 'static, H: Handler<Vec<T>>> Handler<T> for BufferingHandler<T, H> {
    fn handle(&mut self, event: &T) {
        self.buffer.push(event.clone());
        if self.buffer.len() >= self.capacity {
            self.flush();
        }
    }
}

/// Flush any remaining buffered events when the handler is dropped.
/// This fires when the `Collector` exits its run loop and drops all handlers.
impl<T: Clone, H: Handler<Vec<T>>> Drop for BufferingHandler<T, H> {
    fn drop(&mut self) {
        self.flush();
    }
}
