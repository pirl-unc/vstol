use crate::constants::*;
use crate::structs::variant::Variant;
use crate::structs::variant_call::VariantCall;
use crate::structs::variants_list::VariantsList;


/// Verifies that `VariantsList::sort` function sorts variants.
///
/// # Scenario
/// A variants list has variants in the following order: `chr1:10000`, `chr1:100`, and `chr1:1000`.
///
/// # Expected
/// The variants should appear in the order of `chr1:100`, `chr1:1000`, and `chr1:10000`.
///
#[test]
fn sort_variants_list_1() {
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
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
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

    variants_list.add_variant(variant_3);
    variants_list.add_variant(variant_1);
    variants_list.add_variant(variant_2);

    variants_list.sort();

    assert_eq!(variants_list.variants.get(0).unwrap().id, "variant_1");
    assert_eq!(variants_list.variants.get(1).unwrap().id, "variant_2");
    assert_eq!(variants_list.variants.get(2).unwrap().id, "variant_3");
}