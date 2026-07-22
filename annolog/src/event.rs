/// Envelope wrapping user data sent through the collector channel.
///
/// `T` is the user's own type — any struct or enum they define.
/// `Shutdown` is a first-class control signal that stops the collector loop.
pub enum CollectorEvent<T> {
    /// Carries a user-defined value to all registered handlers.
    Data(T),
    /// Signals the collector to drain remaining messages and stop.
    Shutdown,
}
