use crate::constants::*;
use crate::structs::variant::Variant;
use crate::structs::variant_call::VariantCall;
use crate::structs::variants_list::VariantsList;


/// Verifies that `VariantsList::add_variant` function works.
///
/// # Scenario
/// Add two variants to a variants list.
///
/// # Expected
/// The variants list should have 2 variants.
#[test]
fn create_variants_list_1() {
    let mut variants_list: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
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
    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variants_list.add_variant(variant_1);
    variants_list.add_variant(variant_2);
}