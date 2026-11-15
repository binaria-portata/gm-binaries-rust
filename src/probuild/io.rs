//! Output and input primitives emulating the C/C++ runtime of probuild 2.52
//! (re-notes/probuild/core.md section 6).
//!
//! The original is a C++ program (GCC 4.8, static glibc 2.17) that mixes C
//! stdio and C++ iostreams.  Byte-exact output therefore depends on runtime
//! details that this module reproduces:
//!
//! - stdout (`std::cout`, synced with C stdio, fully buffered when it is not
//!   a terminal): collected in a global in-memory buffer (`STDOUT`) and
//!   written at exit (`flush_stdout`).
//! - terminates: where the original ends by SIGSEGV or abort, the port stops
//!   with a clear message and status 1 instead (`fatal`).
//! - stderr (`std::cerr`, C `stderr`): unbuffered, written immediately.  As
//!   stdout is held back, stderr text normally appears before stdout text
//!   when both go to the same terminal or file; this matches the original.
//! - output files (`std::ofstream`, type `OFile`): created or truncated at
//!   open time, content kept in memory and written on `close` or drop.
//!   `exit(1)` (and `fatal`) skips the drops, so a file still open at that
//!   point stays as created (empty).  The original writes through an
//!   8 KiB filebuf, so for files larger than 8 KiB it would already contain
//!   the full blocks; that case is not modelled (core.md section 10).
//! - C++ exceptions: Rust panics carrying `CxxException`, caught in
//!   `catch_main` like main's catch clauses.  Unwinding runs the drops, so
//!   open `OFile`s are flushed, as the C++ destructors flush `ofstream`s.
//! - input: the whole-file loader used for FASTA and `--par` files
//!   (`load_file`) and the two `getline` loop variants (`read_lines`,
//!   `read_lines_eof`).

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::sync::Mutex;

/// Emulated stdout buffer: everything written to `std::cout` so far that has
/// not yet been passed to the real stdout.
static STDOUT: Mutex<Vec<u8>> = Mutex::new(Vec::new());

/// Append to stdout (`std::cout <<`).  Nothing is written to the real stdout
/// until `flush_stdout`, `c_exit` or `fatal`.
pub fn out(b: &[u8]) {
    STDOUT.lock().unwrap().extend_from_slice(b);
}

/// Write to stderr immediately (`std::cerr <<`, `fprintf(stderr, ..)`); both
/// are unbuffered in the original.
pub fn err(b: &[u8]) {
    let mut e = std::io::stderr();
    let _ = e.write_all(b);
}

/// Write the whole emulated stdout buffer to the real stdout and empty it.
/// Called at normal program end, by `c_exit`, and by model modes where the
/// original flushes explicitly (`std::endl` or `flush`), so that the order
/// relative to later stderr output is kept.
pub fn flush_stdout() {
    let mut buffered = STDOUT.lock().unwrap();
    let mut o = std::io::stdout();
    let _ = o.write_all(&buffered);
    let _ = o.flush();
    buffered.clear();
}

/// C `exit(code)` (0x4d2160): stdio buffers are flushed, but C++ objects
/// with automatic storage are not destroyed, so the buffers of open
/// `ofstream`s are lost.  In the port: stdout is written, then the process
/// ends without running the `OFile` drops.
pub fn c_exit(code: i32) -> ! {
    flush_stdout();
    std::process::exit(code)
}

/// A C++ exception thrown inside a handler.  It unwinds (destructors run:
/// open output files are flushed and closed) up to the try/catch of main
/// 0x402ef5, which prints "Out of memory!" for `std::bad_alloc` and "error"
/// for anything else and returns 1 (see `catch_main`).
pub struct CxxException {
    /// True for `std::bad_alloc` (first catch clause, "Out of memory!").
    pub bad_alloc: bool,
}

/// Throw a C++ exception of type `kind` (e.g. `"std::out_of_range"`).  Only
/// `std::bad_alloc` is distinguished by main's catch clauses; `_what` is kept
/// for readability of the call sites and is never printed.
pub fn throw(kind: &str, _what: &str) -> ! {
    std::panic::panic_any(CxxException { bad_alloc: kind == "std::bad_alloc" })
}

/// `std::out_of_range` thrown by `basic_string::substr(pos, n)` when
/// `pos > size()`.  Several modes hit it on unusual definition lines (for
/// example a line consisting of `>` only), which ends the program with
/// "error" and exit code 1.
pub fn throw_substr() -> ! {
    throw("std::out_of_range", "basic_string::substr")
}

/// main's catch clauses (0x40580d `catch (std::bad_alloc&)`, 0x40582a
/// `catch (...)`): run `f`; a `CxxException` panic is turned into the
/// matching stderr message and exit code 1.  Any other panic is a bug in the
/// port and is re-raised unchanged (with the default panic message).
pub fn catch_main<F: FnOnce() -> i32 + std::panic::UnwindSafe>(f: F) -> i32 {
    // Silence the default "thread panicked at" message for emulated C++
    // exceptions only; real Rust panics still print through the old hook.
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        if info.payload().downcast_ref::<CxxException>().is_none() {
            prev(info);
        }
    }));
    match std::panic::catch_unwind(f) {
        Ok(rc) => rc,
        Err(e) => match e.downcast_ref::<CxxException>() {
            Some(x) => {
                err(if x.bad_alloc { b"Out of memory!\n" } else { b"error\n" });
                1
            }
            None => std::panic::resume_unwind(e),
        },
    }
}

/// Stop with a clear error where the original terminates.
///
/// The original ends by SIGSEGV or abort in a few situations (e.g. `vector[0]`
/// on an empty vector, `memset` with a negative length, reads far outside a
/// buffer). The port is deliberately more robust: it writes everything
/// printed so far to stdout, prints `probuild-rs: error: MSG` to stderr and
/// exits with status 1. Open `OFile`s are not flushed (as with `c_exit`), so
/// no half-written output file content is left behind beyond what the
/// original would also have written.
pub fn fatal(msg: &str) -> ! {
    flush_stdout();
    let mut e = std::io::stderr();
    let _ = e.write_all(format!("probuild-rs: error: {msg}\n").as_bytes());
    let _ = e.flush();
    std::process::exit(1)
}

/// `std::ofstream` / `std::fstream` opened for writing.
///
/// The file is created (or truncated, or opened for appending) when the
/// stream is opened, exactly as the original does, so an output file exists
/// even if nothing is ever written to it.  Content is buffered in memory and
/// written on `close` or drop.  Writes to a stream whose open failed are
/// silently ignored, like writes to an `ofstream` in a failed state.
pub struct OFile {
    /// Open file, `None` when not open or when the open failed.
    file: Option<File>,
    /// Content written since opening, not yet passed to the file.
    buf: Vec<u8>,
}

impl OFile {
    /// Not open (like a default-constructed `ofstream`).
    pub fn closed() -> OFile {
        OFile { file: None, buf: Vec::new() }
    }

    /// `ofstream(name, out|trunc)` (`append == false`) or
    /// `ofstream(name, out|app)` (`append == true`).  Check `is_open` for
    /// failure; the original often does not.
    pub fn open(name: &[u8], append: bool) -> OFile {
        let mut o = OFile::closed();
        o.reopen(name, append);
        o
    }

    /// `close()` followed by `open(name, ..)` on the same stream object.
    pub fn reopen(&mut self, name: &[u8], append: bool) {
        self.close();
        let p = std::path::Path::new(std::ffi::OsStr::from_bytes_compat(name));
        let r = if append {
            OpenOptions::new().append(true).create(true).open(p)
        } else {
            OpenOptions::new().write(true).create(true).truncate(true).open(p)
        };
        self.file = r.ok();
    }

    /// `is_open()`.
    pub fn is_open(&self) -> bool {
        self.file.is_some()
    }

    /// `stream << bytes`.  Ignored when the stream is not open.
    pub fn w(&mut self, b: &[u8]) {
        if self.file.is_some() {
            self.buf.extend_from_slice(b);
        }
    }

    /// `stream.flush()`, as done by `std::endl`: pass the buffered content to
    /// the file now. Only observable when the original later terminates (a
    /// termination does not flush), so call it where the original writes `endl`
    /// before a possible abnormal termination.
    pub fn flush(&mut self) {
        if let Some(f) = self.file.as_mut() {
            let _ = f.write_all(&self.buf);
        }
        self.buf.clear();
    }

    /// `close()`: write the buffered content and close the file.  Write
    /// errors are ignored, as the original never checks them.
    pub fn close(&mut self) {
        if let Some(mut f) = self.file.take() {
            let _ = f.write_all(&self.buf);
        }
        self.buf.clear();
    }
}

/// The `ofstream` destructor closes (and thereby flushes) the stream.  This
/// runs at scope end and during unwinding of a `CxxException`, but not for
/// `c_exit` or `fatal`.
impl Drop for OFile {
    fn drop(&mut self) {
        self.close();
    }
}

/// Byte string to `OsStr` (Unix: file names are arbitrary bytes, as in C).
pub trait OsStrBytes {
    fn from_bytes_compat(b: &[u8]) -> &std::ffi::OsStr;
}
impl OsStrBytes for std::ffi::OsStr {
    fn from_bytes_compat(b: &[u8]) -> &std::ffi::OsStr {
        use std::os::unix::ffi::OsStrExt;
        std::ffi::OsStr::from_bytes(b)
    }
}

/// Byte string file name as a `Path`.
pub fn path(b: &[u8]) -> &std::path::Path {
    std::path::Path::new(std::ffi::OsStr::from_bytes_compat(b))
}

/// Whole-file loader 0x40a306 (C stdio): `stat` for the size, `fopen`,
/// `fread` of that many bytes.  Used for every FASTA input and for the
/// `--par` file.  On failure prints one of these messages to stderr and
/// returns `None` (the caller's handler then returns false):
///
/// - empty name: `error: file name is missing`
/// - `stat` fails: `error: file 'NAME' is not found`
/// - `fopen` fails (e.g. no permission): `error: on opening file 'NAME'`
/// - read error (e.g. NAME is a directory): `error: on reading file 'NAME'`
///
/// At most the size reported by `stat` is read; a shorter read (file shrank)
/// keeps what was read.
pub fn load_file(name: &[u8]) -> Option<Vec<u8>> {
    if name.is_empty() {
        err(b"error: file name is missing\n");
        return None;
    }
    let report = |prefix: &str, suffix: &str| {
        let mut m = prefix.as_bytes().to_vec();
        m.extend_from_slice(name);
        m.extend_from_slice(suffix.as_bytes());
        err(&m);
    };
    let p = path(name);
    let md = match std::fs::metadata(p) {
        Ok(m) => m,
        Err(_) => {
            report("error: file '", "' is not found\n");
            return None;
        }
    };
    let mut f = match File::open(p) {
        Ok(f) => f,
        Err(_) => {
            report("error: on opening file '", "'\n");
            return None;
        }
    };
    let size = md.len() as usize;
    let mut v = vec![0u8; size];
    let mut got = 0;
    use std::io::Read;
    while got < size {
        match f.read(&mut v[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(_) => {
                report("error: on reading file '", "'\n");
                return None;
            }
        }
    }
    v.truncate(got);
    Some(v)
}

/// The loop `while (std::getline(f, line)) { .. }` over a file opened as
/// `std::ifstream` (sid tables, `--fasta_list`).  Returns `None` if the open
/// fails.  Lines are returned without their '\n'; a final empty piece after
/// a trailing newline is NOT included (that `getline` call fails), but a
/// last line without a newline is.  '\r' is kept as part of the line.
pub fn read_lines(name: &[u8]) -> Option<Vec<Vec<u8>>> {
    Some(split_lines(&read_ifstream(name)?))
}

/// The loop `while (!f.eof()) { std::getline(f, line); .. }` (GFF regions of
/// `--cut`, the `--mask` file): every '\n'-separated piece, INCLUDING the
/// (possibly empty) piece after the last '\n', because `eof()` only becomes
/// true after the `getline` that hits the end.  An empty file gives one empty
/// line.  The callers skip or ignore empty lines, but the difference matters
/// for their error messages.
pub fn read_lines_eof(name: &[u8]) -> Option<Vec<Vec<u8>>> {
    let data = read_ifstream(name)?;
    Some(data.split(|&c| c == b'\n').map(|s| s.to_vec()).collect())
}

/// Content of a file opened with `std::ifstream`.  A directory opens
/// successfully but reading fails, which behaves like an empty file.
fn read_ifstream(name: &[u8]) -> Option<Vec<u8>> {
    use std::io::Read;
    let mut f = File::open(path(name)).ok()?;
    let mut data = Vec::new();
    if f.read_to_end(&mut data).is_err() {
        data.clear();
    }
    Some(data)
}

/// Split `data` into the lines `while (getline(..))` delivers: pieces ending
/// in '\n' (without it), plus a non-empty last piece without '\n'.
pub fn split_lines(data: &[u8]) -> Vec<Vec<u8>> {
    let mut lines = Vec::new();
    let mut start = 0;
    for (i, &c) in data.iter().enumerate() {
        if c == b'\n' {
            lines.push(data[start..i].to_vec());
            start = i + 1;
        }
    }
    if start < data.len() {
        lines.push(data[start..].to_vec());
    }
    lines
}
