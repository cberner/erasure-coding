use crate::gf256::Polynomial;

pub struct Encoder {
    data_blocks: u8,
    repair_blocks: u8
}

impl Encoder {
    pub fn new(data_blocks: u8, repair_blocks: u8) -> Encoder {
        Encoder {
            data_blocks,
            repair_blocks
        }
    }

    // Returns tuple of (data blocks, repair blocks)
    pub fn encode(&self, data: &[u8]) -> (Vec<u8>, Vec<u8>) {
        // TODO: support blocks with more than 1 byte
        assert_eq!(data.len(), self.data_blocks as usize);
        let mut coefficients = data.to_vec();
        coefficients.extend(vec![0; self.repair_blocks as usize]);

        let poly = Polynomial::new(&coefficients);
        let generator_polynomial = Polynomial::create_generator_polynomial(self.repair_blocks);
        let (_, remainder) = poly.div(&generator_polynomial);

        return (data.to_vec(), remainder.into_coefficients())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn encode() {
        let data = vec![1, 2, 3];
        let expected_repair = vec![4, 4];

        let encoder = super::Encoder::new(data.len() as u8, expected_repair.len() as u8);
        let (encoded_data, repair) = encoder.encode(&data);

        assert_eq!(data, encoded_data);
        assert_eq!(expected_repair, repair);
    }
}
