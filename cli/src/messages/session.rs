use serde_json::Value as JsonValue;

pub fn subject(session_name: &str) -> String {
    format!("Session: {}", session_name)
}

pub fn invite_body(owner: &str, session_name: &str, session_id: &str) -> String {
    format!(
        "{} invited you to the session \"{}\" (ID: {}).",
        owner, session_name, session_id
    )
}

pub fn invite_metadata(
    session_id: &str,
    session_name: &str,
    owner: &str,
    description: &Option<String>,
    created_at_rfc3339: &str,
) -> JsonValue {
    serde_json::json!({
        "session_invite": {
            "session_id": session_id,
            "session_name": session_name,
            "from": owner,
            "description": description,
            "created_at": created_at_rfc3339,
        },
        // Also include session_chat metadata so the invite anchors the session thread.
        "session_chat": {
            "session_id": session_id,
            "session_name": session_name,
            "from": owner,
            "created_at": created_at_rfc3339,
        }
    })
}

pub fn invite_response_body(
    responder: &str,
    session_name: &str,
    session_id: &str,
    accepted: bool,
    reason: &Option<String>,
) -> String {
    if accepted {
        format!(
            "{} accepted your session invite for \"{}\" (ID: {}).",
            responder, session_name, session_id
        )
    } else {
        format!(
            "{} declined your session invite for \"{}\" (ID: {}).{}",
            responder,
            session_name,
            session_id,
            reason
                .as_ref()
                .map(|r| format!(" Reason: {}", r))
                .unwrap_or_default()
        )
    }
}

pub fn invite_response_metadata(
    session_id: &str,
    session_name: &str,
    responder: &str,
    accepted: bool,
    reason: &Option<String>,
    responded_at_rfc3339: &str,
    created_at_rfc3339: &str,
) -> JsonValue {
    let status = if accepted { "accepted" } else { "rejected" };
    serde_json::json!({
        "session_invite_response": {
            "session_id": session_id,
            "status": status,
            "responder": responder,
            "reason": reason,
            "session_name": session_name,
            "responded_at": responded_at_rfc3339,
        },
        "session_chat": {
            "session_id": session_id,
            "session_name": session_name,
            "from": responder,
            "created_at": created_at_rfc3339,
        }
    })
}

pub fn chat_metadata(
    session_id: &str,
    session_name: &str,
    from: &str,
    created_at_rfc3339: &str,
) -> JsonValue {
    serde_json::json!({
        "session_chat": {
            "session_id": session_id,
            "session_name": session_name,
            "from": from,
            "created_at": created_at_rfc3339,
        }
    })
}
