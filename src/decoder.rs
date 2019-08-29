use crate::gf256::{Polynomial, Octet, mulassign_scalar};
use crate::base::Block;
use crate::block_polynomial::BlockPolynomial;

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

    fn calculate_syndrome(&self, data: &[Option<Block>], repair: &[Block]) -> BlockPolynomial {
        let block_length = repair[0].len();
        let mut full_with_zeros: Vec<Block> = data.iter().map(|x| x.clone().unwrap_or(vec![0; block_length])).collect();
        full_with_zeros.extend(repair.to_vec());
        let full_poly = BlockPolynomial::new(full_with_zeros);
        // Pad with a zero
        let mut syndrome = vec![vec![0; block_length]];
        for i in 0..self.repair_blocks {
            syndrome.push(full_poly.eval(&Octet::alpha(i)));
        }
        syndrome.reverse();
        return BlockPolynomial::new(syndrome);
    }

    fn calculate_erasure_locator(&self, erasures: &[bool]) -> Polynomial {
        let mut locator = Polynomial::new(&[1]);
        for i in 0..erasures.len() {
            if erasures[i] {
                let location = self.data_blocks + self.repair_blocks - 1 - i as u8;
                locator = locator.mul(&Polynomial::new(&[Octet::alpha(location as u8).byte(), 1]));
            }
        }
        return locator;
    }

    fn calculate_erasure_evaluator(&self, syndrome: &BlockPolynomial, erasure_locator: &Polynomial) -> BlockPolynomial {
        let mut poly = Polynomial::new(&vec![0; (erasure_locator.coefficients.len() + 1) as usize]);
        poly.coefficients[0] = Octet::one();

        return syndrome.mul_poly(erasure_locator).zero_extend_div_remainder(0, &poly);
    }

    // Forney algorithm
    fn calculate_delta_correction(&self, erasure_evaluator: &BlockPolynomial, erasures: &[bool]) -> BlockPolynomial {
        let block_length = erasure_evaluator.coefficient_arrays[0].len();
        let mut error_locations = vec![];
        let mut full_message_erasure_positions = vec![];
        for i in 0..erasures.len() {
            if erasures[i] {
                full_message_erasure_positions.push(i);
                let location = self.data_blocks + self.repair_blocks - 1 - i as u8;
                let l = 255 - location as i32;
                error_locations.push(Octet::new(2).pow(-l));
            }
        }

        let mut error_magnitudes = vec![vec![0; block_length]; self.data_blocks as usize];
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

            let mut y = erasure_evaluator.eval(&inverse_error_value);
            mulassign_scalar(&mut y, &(error_value / &locator_prime));

            error_magnitudes[full_message_erasure_positions[i]] = y;
        }

        return BlockPolynomial::new(error_magnitudes);
    }

    // TODO: support missing repair blocks
    // Returns tuple of (data blocks, repair blocks)
    pub fn decode(&self, data: &[Option<Block>], repair: &[Block]) -> Vec<u8> {
        assert_eq!(data.len(), self.data_blocks as usize);
        assert_eq!(repair.len(), self.repair_blocks as usize);
        let block_length = repair[0].len();
        let data_with_zeros: Vec<Block> = data.iter().map(|x| x.clone().unwrap_or(vec![0; block_length])).collect();

        let syndrome = self.calculate_syndrome(data, repair);
        let erasures: Vec<bool> = data.iter().map(|x| x.is_none()).collect();
        let locator = self.calculate_erasure_locator(&erasures);
        let erasure_evaluator = self.calculate_erasure_evaluator(&syndrome, &locator);
        let correction = self.calculate_delta_correction(&erasure_evaluator, &erasures);

        let mut data_with_zeros_poly = BlockPolynomial::new(data_with_zeros);
        data_with_zeros_poly.addassign(&correction);

        let repaired_data = data_with_zeros_poly.into_blocks();
        let mut repaired = Vec::with_capacity(self.data_blocks as usize * block_length);
        for j in 0..self.data_blocks as usize {
            repaired.extend(repaired_data[j].clone());
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
        let erasures = vec![true, false, false, true, true];
        assert_eq!(Polynomial::new(&[19, 28, 152, 1]),
                   Decoder::new(5, 3).calculate_erasure_locator(&erasures));
    }
}
