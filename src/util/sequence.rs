
use simple_error::{bail, SimpleError};

/// Creates a reverse complement sequence from an input.
/// # Arguments
/// * `original` - the sequence to rev-comp
/// # Errors
/// * if any non-ACGNT character is provided
pub fn reverse_complement(original: &[u8]) -> Result<Vec<u8>, SimpleError> {
    original.iter()
        .rev()
        .map(|c| {
            match c {
                b'A' => Ok(b'T'),
                b'C' => Ok(b'G'),
                b'G' => Ok(b'C'),
                b'T' => Ok(b'A'),
                b'N' => Ok(b'N'),
                _ => bail!("Unexpected character for reverse-complement: {c}")
            }
        })
        .collect::<Result<Vec<u8>, SimpleError>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let sequence = b"ACCGGGTN";
        let expected = b"NACCCGGT";
        let rc_result = reverse_complement(sequence).unwrap();
        assert_eq!(&rc_result, expected);
    }

    #[test]
    fn test_reverse_complement_invalid() {
        let sequence = b"b";
        let rc_result = reverse_complement(sequence);
        assert!(rc_result.is_err());
    }
}