#![windows_subsystem = "windows"]
use fltk::{app::*, button::*, dialog::*, input::*, output::*, prelude::*, window::*};
use statrs::function::gamma::*;

#[derive(Debug, Clone, Copy)]
pub enum Message {
    Add,
    Clear,
}

#[derive(Clone, Debug)]
// Define a struct for the form fields
struct Controls {
    n_in: IntInput,
    es_in: FloatInput,
    p_in: FloatInput,
    n_out: Output,
    es_out: Output,
    p_out: Output,
}

#[derive(Clone, Debug)]
// Define a struct for the form fields
struct Data {
    n: i32,
    n_data: Vec<i32>,
    es_data: Vec<f64>,
    p_data: Vec<f64>,
}

fn main() {
    // Get app handle
    let app = App::default();

    // Main Window
    let mut wind = Window::new(
        100,
        100,
        300,
        300,
        "Meta-Analysis Calculator v1.0",
    );

    // Define Buttons
    let mut add_button = Button::new(30, 240, 100, 40, "Add Trial");
    let mut clear_button = Button::new(170, 240, 100, 40, "Clear");

    // Define Text Controls
    let mut controls = Controls {
        n_in: IntInput::new(100, 30, 75, 22, "Sample Size"),
        es_in: FloatInput::new(100, 55, 75, 22, "Effect Size"),
        p_in: FloatInput::new(100, 80, 75, 22, "p-Value"),
        n_out: Output::new(100, 135, 75, 22, "# Trials"),
        es_out: Output::new(100, 160, 75, 22, "Effect Size"),
        p_out: Output::new(100, 185, 75, 22, "p-Value"),
    };

    // Make a place to store data
    let mut data = Data {
        n: 0,
        n_data: Vec::new(),
        es_data: Vec::new(),
        p_data: Vec::new(),
    };

    // Show the window
    wind.end();
    wind.show();

    // Setup the message handler
    let (s, r) = channel::<Message>();

    // Attach messages to event emitters
    add_button.emit(s, Message::Add);
    clear_button.emit(s, Message::Clear);

    // Main Message Loop
    while app.wait() {
        if let Some(msg) = r.recv() {
            match msg {
                Message::Add => add(&mut controls, &mut data),
                Message::Clear => clear(&mut controls, &mut data),
            }
        }
    }
}

// Handle add button
fn add(controls: &mut Controls, data: &mut Data) {
    // Parse and push data into vectors
    let n = match controls.n_in.value().parse::<i32>() {
        Ok(v) => v,
        Err(_) => return,
    };
    let e = match controls.es_in.value().parse::<f64>() {
        Ok(v) => v,
        Err(_) => return,
    };
    let p = match controls.p_in.value().parse::<f64>() {
        Ok(v) => v,
        Err(_) => return,
    };

    // Cannot have less than 1 sample count
    if n < 1 {
        alert(368, 265,"Invalid Sample Size, must be > 0");
        return;
    }

    // p value must be 0-1
    if p <= 0.0 || p >= 1.0 {
        alert(368, 265,"Invalid P-Value, must be 0 > p < 1");
        return;
    }

    // Keep the validated data
    data.n_data.push(n);
    data.es_data.push(e);
    data.p_data.push(p);

    // Clear the input fields
    controls.n_in.set_value(&"");
    controls.es_in.set_value(&"");
    controls.p_in.set_value(&"");

    // Increment the number of trials count
    data.n += 1;

    // Display data
    controls.n_out.set_value(&format!("{}", data.n));
    controls.es_out.set_value(&science_pretty_format(es(&data.es_data, &data.n_data), 3));

    match p_value(&data.p_data) {
        Ok(v) => controls.p_out.set_value(&science_pretty_format(v, 3)),
        Err(e) => {
            alert(368, 265, &e);
        }
    }
}

// Handle clear button
fn clear(controls: &mut Controls, data: &mut Data) {
    //Clear Data
    data.n = 0;
    data.n_data.clear();
    data.es_data.clear();
    data.p_data.clear();

    //Clear form fields
    controls.n_in.set_value(&"");
    controls.es_in.set_value(&"");
    controls.p_in.set_value(&"");
    controls.n_out.set_value(&"");
    controls.es_out.set_value(&"");
    controls.p_out.set_value(&"");
}

// Calculate meta p-value
fn p_value(p_data: &[f64]) -> Result<f64, String> {
    let mut sum: f64 = 0.0;

    for v in p_data.iter() {
        sum += -2.0 * v.ln();
    }

    chisqr(p_data.len() * 2,sum)
}

// Calculate p-Value from chi-square sum
fn chisqr(dof: usize, cv: f64) -> Result<f64, String> {
    if cv < 0.0 {
        return Err("Invalid P-Value, must be 0 > p < 1".to_string());
    }

    if dof < 1 {
        return Err("DOF must be 1 or greater".to_string());
    }

    let k: f64 = dof as f64 * 0.5;
    let x = cv * 0.5;

    if dof == 2 {
        return Ok((-1.0 * x).exp());
    }

    Ok(1.0 - (log_igf(k, x) - ln_gamma(k)).exp())
}

// Calculate ES value
fn es(es: &[f64], n: &[i32]) -> f64 {
    let mut numerator: f64 = 0.0;
    let mut denominator: f64 = 0.0;

    for (i, trial) in es.iter().enumerate() {
        numerator += n[i] as f64 * trial;
        denominator += n[i] as f64;
    }

    numerator / denominator
}

// Calculate KM
fn km(s: f64, z: f64) -> f64 {
    let mut sum = 1.0;
    let mut nom = 1.0;
    let mut denom = 1.0;
    let mut s = s;
    let mut tsum: f64;

    loop {
        nom *= z;
        s += 1.0;
        denom *= s;
        tsum = sum + nom / denom;
        if tsum - sum <= f64::EPSILON {
            break;
        }
        sum = tsum;
    }

    sum
}

// Calculate log of incomplete gamma function
fn log_igf(s: f64, z: f64) -> f64 {
    if z < 0.0 {
        return 0.0;
    }

    let sc = z.ln() * s - z - s.ln();
    let k = km(s, z);

    k.ln() + sc
}

// Pretty Format Scientific Numbers
fn science_pretty_format(value: f64, digits: usize) -> String {
    if value.abs() == 0.0 {
        "0".to_string();
    }
    if value.abs() >= 10000.0 || value.abs() < 0.001 {
        format!("{:.*e}", digits, value);
    }
    format!("{:.*}", digits, value)
        .trim_end_matches(|c| c == '0')
        .trim_end_matches(|c| c == '.')
        .to_string()
}
