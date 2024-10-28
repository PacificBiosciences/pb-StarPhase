
/// This is an amalgamation of the Multinomial distribution from statrs and the factorial function which needs to be in ln mode.
/// Given a set of expected probabilities and observations, it will return the log-likelihood of that observation distribution.
/// References: https://github.com/statrs-dev/statrs/blob/e8e9c61b860241c70f9c71f2e07fbd6dde2cf44f/src/function/factorial.rs#L71
///             https://github.com/statrs-dev/statrs/blob/e8e9c61b860241c70f9c71f2e07fbd6dde2cf44f/src/distribution/multinomial.rs#L222
/// # Arguments
/// * `probs` - The probability of observing each category, the sum is expected to add to 1.0
/// * `obs` - The number of times each category was observed
/// # Panics
/// * if `probs.len() != obs.len()`
pub fn multinomial_ln_pmf(probs: &[f64], obs: &[u64]) -> f64 {
    use statrs::function::factorial::ln_factorial;

    // copy of sanity checks
    if probs.len() != obs.len() {
        panic!("Expected probs and obs to have equal lengths.");
    }

    // we just derive total count here
    let total_count: u64 = obs.iter().sum();
    assert!(total_count > 0);

    // here is where have to implement the multinomial in log space
    let mut coeff = ln_factorial(total_count); // initialized to total_count!
    for observed_count in obs.iter() {
        // subtract out each observed_count! (this is division in normal space)
        coeff -= ln_factorial(*observed_count);
    }

    // the factorial just gets added because it is multiplied in combinatorial space
    let val = coeff
        + probs.iter() // this whole zip-iter is just copy pasted, but it's basically doing all the probability multiplication
            .zip(obs.iter())
            .map(|(pi, xi)| *xi as f64 * pi.ln())
            .fold(0.0, |acc, x| acc + x);
    val
}

#[cfg(test)]
mod tests {
    use super::*;

    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_multinomial() {
        // 1 category - should basically always be 1.0
        let probs = [1.0];
        let obs = [10];
        assert_approx_eq!(multinomial_ln_pmf(&probs, &obs), 1.0_f64.ln());

        // 2 categories (binomial)
        let probs = [0.25, 0.75];
        let obs = [1, 3]; // normal ratios
        let expected_prob: f64 = 4.0 * 0.25 * 0.75_f64.powf(3.0);
        assert_approx_eq!(multinomial_ln_pmf(&probs, &obs), expected_prob.ln());

        let obs = [3, 1]; // abnormal ratios
        let expected_prob = 4.0 * 0.25_f64.powf(3.0) * 0.75;
        assert_approx_eq!(multinomial_ln_pmf(&probs, &obs), expected_prob.ln());

        // 3+ categories
        let probs = [0.25, 0.25, 0.5];
        let obs = [1, 1, 2]; // normal ratios
        let expected_prob: f64 = (4.0 * 3.0 * 2.0 / 2.0) * 0.25 * 0.25 * 0.5_f64.powf(2.0);
        assert_approx_eq!(multinomial_ln_pmf(&probs, &obs), expected_prob.ln());

        let obs = [2, 2, 0]; // abnormal ratios
        let expected_prob: f64 = (4.0 * 3.0 * 2.0 / 2.0 / 2.0) * 0.25_f64.powf(4.0);
        assert_approx_eq!(multinomial_ln_pmf(&probs, &obs), expected_prob.ln());
    }
}