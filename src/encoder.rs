use crate::gf256::Polynomial;
use crate::base::Block;

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
    pub fn encode(&self, data: &[u8]) -> (Vec<Block>, Vec<Block>) {
        assert_eq!(data.len() % self.data_blocks as usize, 0);
        let block_length = data.len() / self.data_blocks as usize;
        // Data stripped across blocks: [0, block_length) in the first block,
        // then [block_length, 2 * block_length) in the second
        let mut data_blocks = vec![];
        for i in 0..self.data_blocks as usize {
            data_blocks.push(data[i * block_length..(i + 1)*block_length].to_vec())
        }

        let mut repair_blocks = vec![Vec::with_capacity(block_length); self.repair_blocks as usize];
        let generator_polynomial = Polynomial::create_generator_polynomial(self.repair_blocks);

        // Allocate this outside the loop to avoid excess memory allocations
        let mut coefficients = vec![0; (self.data_blocks + self.repair_blocks) as usize];
        for i in 0..block_length {
            // Encode one byte from each block, so that the polynomial spans all the blocks
            for j in 0..self.data_blocks as usize {
                coefficients[j] = data_blocks[j][i];
            }
            // Set all repair coefficients to zero
            for j in 0..self.repair_blocks as usize {
                coefficients[self.data_blocks as usize + j] = 0;
            }

            let poly = Polynomial::new(&coefficients);
            let (_, remainder) = poly.div(&generator_polynomial);
            for j in 0..self.repair_blocks as usize {
                repair_blocks[j].push(remainder.coefficients[j].byte());
            }
        }

        return (data_blocks, repair_blocks);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn encode() {
        let data = vec![1, 2, 3];
        let expected_data = vec![vec![1], vec![2], vec![3]];
        let expected_repair = vec![vec![4], vec![4]];

        let encoder = super::Encoder::new(data.len() as u8, expected_repair.len() as u8);
        let (encoded_data, repair) = encoder.encode(&data);

        assert_eq!(expected_data, encoded_data);
        assert_eq!(expected_repair, repair);
    }
}
