use crate::gf256::{Polynomial, Octet};
use crate::base::Block;

pub struct Decoder {
    data_blocks: u8,
    repair_blocks: u8
}

impl Decoder {
    pub fn new(data_blocks: u8, repair_blocks: u8) -> Decoder {
        Decoder {
            data_blocks,
            repair_blocks
        }
    }

    fn calculate_syndrome(&self, data: &[Option<u8>], repair: &[u8]) -> Polynomial {
        let mut full_with_zeros: Vec<u8> = data.iter().map(|x| x.unwrap_or(0)).collect();
        full_with_zeros.extend(repair);
        let full_poly = Polynomial::new(&full_with_zeros);
        // Pad with a zero
        let mut syndrome = vec![0; self.repair_blocks as usize + 1];
        for i in 0..self.repair_blocks {
            syndrome[i as usize + 1] = full_poly.eval(&Octet::alpha(i)).byte();
        }
        syndrome.reverse();
        return Polynomial::new(&syndrome);
    }

    fn calculate_erasure_locator(&self, data: &[Option<u8>]) -> Polynomial {
        let mut locator = Polynomial::new(&[1]);
        for i in 0..data.len() {
            if data[i].is_none() {
                let location = self.data_blocks + self.repair_blocks - 1 - i as u8;
                locator = locator.mul(&Polynomial::new(&[Octet::alpha(location as u8).byte(), 1]));
            }
        }
        return locator;
    }

    fn calculate_erasure_evaluator(&self, syndrome: &Polynomial, erasure_locator: &Polynomial) -> Polynomial {
        let mut poly = Polynomial::new(&vec![0; (erasure_locator.coefficients.len() + 1) as usize]);
        poly.coefficients[0] = Octet::one();

        let (_, remainder) = syndrome.mul(erasure_locator).div(&poly);
        return remainder
    }

    // Forney algorithm
    fn correct_erasures(&self, syndrome: &Polynomial, data: &[Option<u8>], repair: &[u8]) -> Polynomial {
        let locator = self.calculate_erasure_locator(data);
        let erasure_evaluator = self.calculate_erasure_evaluator(&syndrome, &locator);

        let mut error_locations = vec![];
        let mut full_message_erasure_positions = vec![];
        for i in 0..data.len() {
            if data[i].is_none() {
                full_message_erasure_positions.push(i);
                let location = self.data_blocks + self.repair_blocks - 1 - i as u8;
                let l = 255 - location as i32;
                error_locations.push(Octet::new(2).pow(-l));
            }
        }

        let mut error_magnitudes = vec![Octet::zero(); (self.data_blocks + self.repair_blocks) as usize];
        let error_count = error_locations.len();
        for (i, error_value) in error_locations.iter().enumerate() {
            let inverse_error_value = &Octet::one() / error_value;

            let mut locator_prime = Octet::one();
            for j in 0..error_count {
                if j != i {
                    locator_prime = &locator_prime * &(Octet::one() - &inverse_error_value * &error_locations[j]);
                }
            }
            assert_ne!(locator_prime, Octet::zero());

            let y = error_value * &erasure_evaluator.eval(&inverse_error_value);
            let magnitude = y / locator_prime;

            error_magnitudes[full_message_erasure_positions[i]] = magnitude;
        }

        let mut full_with_zeros: Vec<u8> = data.iter().map(|x| x.unwrap_or(0)).collect();
        full_with_zeros.extend(repair);
        let full_poly = Polynomial::new(&full_with_zeros);
        let error_bytes: Vec<u8> = error_magnitudes.iter().map(Octet::byte).collect();
        return full_poly.add(&Polynomial::new(&error_bytes));
    }

    // TODO: support missing repair blocks
    // Returns tuple of (data blocks, repair blocks)
    pub fn decode(&self, data: &[Option<Block>], repair: &[Block]) -> Vec<u8> {
        assert_eq!(data.len(), self.data_blocks as usize);
        assert_eq!(repair.len(), self.repair_blocks as usize);
        let block_length = repair[0].len();
        let mut repaired = vec![0; self.data_blocks as usize * block_length];

        // Allocate this outside the loop to avoid excess memory allocations
        let mut data_coefficients = vec![None; self.data_blocks as usize];
        let mut repair_coefficients = vec![0; self.repair_blocks as usize];
        for i in 0..block_length {
            // Take the i'th byte out of each block, since the polynomials span the blocks
            for j in 0..self.data_blocks as usize {
                if let Some(ref data_block) = data[j] {
                    data_coefficients[j] = Some(data_block[i]);
                } else {
                    data_coefficients[j] = None;
                }
            }
            for j in 0..self.repair_blocks as usize {
                repair_coefficients[j] = repair[j][i];
            }
            let syndrome = self.calculate_syndrome(&data_coefficients, &repair_coefficients);
            let poly = self.correct_erasures(&syndrome, &data_coefficients, &repair_coefficients);
            let repaired_data = poly.into_coefficients();
            for j in 0..self.data_blocks as usize {
                repaired[j*block_length + i] = repaired_data[j];
            }
        }

        return repaired;
    }
}

#[cfg(test)]
mod tests {
    use crate::Decoder;
    use crate::gf256::Polynomial;

    #[test]
    fn erasure_locator() {
        let data = vec![None, Some(1), Some(2), None, None];
        assert_eq!(Polynomial::new(&[19, 28, 152, 1]),
                   Decoder::new(5, 3).calculate_erasure_locator(&data));
    }
}
