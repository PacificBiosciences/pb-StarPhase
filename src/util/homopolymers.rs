
use itertools::Itertools;

/// This will take a sequence and homo-polymer compress it, remove many common errors but also potentially HPC true variation.
/// # Arguments
/// * `sequence` - the sequence to run HPC on
/// # Errors
/// * if String::from_utf8 fails
pub fn hpc(sequence: &str) -> Result<String, std::string::FromUtf8Error> {
    String::from_utf8(
        hpc_bytes(sequence.as_bytes())
    )
}

/// This is a generic HPC for any type that can be compared & cloned
/// # Arguments
/// * `sequence` - the sequence to run HPC on
pub fn hpc_bytes<T: std::cmp::PartialEq + Clone>(sequence: &[T]) -> Vec<T> {
    sequence.iter()
        .dedup()
        .cloned()
        .collect::<Vec<T>>()
}

/// This will take a sequence and an index and return the new index in a homo-polymer compressed version.
/// # Arguments
/// * `sequence` - the sequence to run HPC on
/// * `position` - the position we want to find in the new HPC
pub fn hpc_pos<T: std::cmp::PartialEq>(sequence: &[T], position: usize) -> usize {
    let mut total_length = 0;
    let mut offset = 0;

    for (l, _v) in sequence.iter().dedup_with_count() {
        total_length += l;
        if position < total_length {
            break;
        }
        offset += 1;
    }

    offset
}


/// This will take a sequence and homo-polymer compress it, remove many common errors but also potentially HPC true variation.
/// It then aligns it to an HPC guide sequence and adds a prefix buffer if needed.
/// # Arguments
/// * `sequence` - the sequence to run HPC on
/// * `guide_sequence` - the sequence we want to then align against and add a prefix for the offset, this get HPC before we do the alignment
/// * `guide_offset` - the offset into the guide where the sequence starts to align
/// # Errors
/// * if String::from_utf8 fails for either sequence
pub fn hpc_with_guide(sequence: &str, guide_sequence: &str, guide_offset: usize) -> Result<(String, usize), Box<dyn std::error::Error>> {
    // first, HPC the input
    let hpc_sequence = hpc(sequence)?;

    // now figure out where this should start in the guide sequence
    let hpc_guide_offset = hpc_pos(guide_sequence.as_bytes(), guide_offset);
    
    // now solve the prefix and append it
    // let prefix = String::from_utf8(vec![b'*'; hpc_guide_offset])?;
    // let final_hpc = prefix + &hpc_sequence;
    Ok((hpc_sequence, hpc_guide_offset))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hpc() {
        let sequence = "AACAAAAAAGGGTAACAA";
        let expected = "ACAGTACA";
        let hpc_result = hpc(sequence).unwrap();
        assert_eq!(&hpc_result, expected);
    }

    #[test]
    fn test_hpc_pos() {
        let sequence = "AACCCGTTTT";
        for (i, &c) in sequence.as_bytes().iter().enumerate() {
            let expected = match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("should not happen")
            };
            assert_eq!(expected, hpc_pos(sequence.as_bytes(), i));
        }
    }

    #[test]
    fn test_hpc_guide() {
        let guide_sequence = "ATTGGGGGAACCCGTTTT";
        let sequence =              "GAACCCGTTTT";
        let hpc_result = hpc_with_guide(sequence, guide_sequence, 6).unwrap();
        assert_eq!(hpc_result.0, "GACGT"); // ATT -> **
        assert_eq!(hpc_result.1, 2);
    }
}