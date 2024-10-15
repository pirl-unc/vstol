use std::collections::HashMap;
use crate::constants::*;
use crate::structs::genomic_range::GenomicRange;
use crate::structs::genomic_ranges_list::GenomicRangesList;
use crate::structs::variant::Variant;
use crate::structs::variant_call::VariantCall;
use crate::structs::variants_list::VariantsList;


/// Verifies that `VariantsList::overlap` function correctly identifies variants overlapping a
/// GenomicRangeList object.
///
/// # Scenario
/// A variants list has 3 variants but only one of them overlaps the input GenomicRangesList object.
///
/// # Expected
/// 1 variant call in the overlapping variant call IDs list.
#[test]
fn overlap_variants_list_1() {
    let mut variants_list: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        INSERTION.to_string(),
        "T".to_string(),
        "TAAA".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        3,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.position_1 = 1000;
    variant_call_2.position_2 = 1000;
    let mut variant_call_3: VariantCall = variant_call_1.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_3.position_1 = 10000;
    variant_call_3.position_2 = 10000;
    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list.add_variant(variant_1);
    variants_list.add_variant(variant_2);
    variants_list.add_variant(variant_3);

    let mut granges_list: GenomicRangesList = GenomicRangesList::new();
    let grange_1: GenomicRange = GenomicRange::new(
        "chr1",
        1,
        200
    );
    granges_list.add_genomic_range(&grange_1);

    let overlapping_variant_call_ids: HashMap<String,Vec<String>> = variants_list.overlap(
        granges_list,
        1
    );

    assert_eq!(overlapping_variant_call_ids.keys().len(), 1);
    assert_eq!(overlapping_variant_call_ids.contains_key("variant_call_1"), true);
    assert_eq!(overlapping_variant_call_ids.contains_key("variant_call_2"), false);
    assert_eq!(overlapping_variant_call_ids.contains_key("variant_call_3"), false);
}
